
#获得各个组织的增强子的bed4文件，为igv输入。
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if(NR>1&&$NF==tissue){print $1,$2,$3,$4} }' enh_all_count_spe_melt >/media/ding/000B49000006264C/eRNA_project/figure/table/for_igv/$i"_enh.bed"
  echo $i"完成"
done

#获得各个组织的IGV要提取的范围：上下游10000bp:
#IGV session的xml文件先建立一个组织的，然后替换组织名，批量修改。

tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  #设置要使用igv提取的范围：上下游10kb:
  cd /home/ding/all_cmd/script/enh_statistics
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if(NR>1&&$NF==tissue){print $1,$2-10000,$3+10000,$4} }' enh_all_count_spe_melt >/media/ding/000B49000006264C/eRNA_project/figure/table/for_igv/$i"_enh_interval.bed"
  
  #使用bedToIgv获得每个增强子附近的修饰特征。
  cd  /media/ding/000B49000006264C/eRNA_project/figure/table/for_igv
  bedToIgv -i /media/ding/000B49000006264C/eRNA_project/figure/table/for_igv/$i"_enh_interval.bed"   -sess /media/ding/000B49000006264C/eRNA_project/figure/table/for_igv/$i".xml" -path /media/ding/000B49000006264C/eRNA_project/figure/table/for_igv/$i -name  >/media/ding/000B49000006264C/eRNA_project/figure/table/for_igv/$i".igv.batch"
  
  #对于每个batch的图片名称进行修改，去掉位置，只保留enh_1.png类似的文件名。
  cd  /media/ding/000B49000006264C/eRNA_project/figure/table/for_igv
  awk 'BEGIN{OFS=" "}{ if(NR<=2){print } if(NR>2&&NR%2!=0){print } if(NR>2&&NR%2==0){split($2,a,"_");print $1,a[4]"_"a[5]} }' /media/ding/000B49000006264C/eRNA_project/figure/table/for_igv/$i".igv.batch">./a
  cp ./a /media/ding/000B49000006264C/eRNA_project/figure/table/for_igv/$i".igv.batch"
  echo $i"完成"
  rm ./a
done


#再次注意：先设置的是Adrenal的igv session:
#先设置好Adrenal的IGV session的xml文件，然后批量生成各个组织的xml文件，然后替换组织。
cd  /media/ding/000B49000006264C/eRNA_project/figure/table/for_igv
tissue=(Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  #全部拷贝： 
  cp ./Adrenal.xml ./$i".xml"
  #各自替换
  sed -i 's/Adrenal/'"$i"'/g' ./$i".xml"
  echo $i"完成" 
done


#生成组织特异性的调控数据：作为数据库ts搜素原始文件：
#数据来源:TS_enh_TF_targetgene_3  没有将TF和gene同名的基因修改为gene_
awk 'BEGIN{OFS="\t";print "enh_name","isbid","TF","tissue","gene_name","TF_exp","gene_TS_score","TF_TS_score"}{ print $0 }' /home/ding/all_cmd/script/enh_statistics/TS_enh_TF_targetgene_3 > /media/ding/000B49000006264C/eRNA_project/figure/table/database/TS_enh_TF_targetgene_3


#修改生成特异性网络的联动菜单.js文件
#awk '' /home/ding/all_cmd/script/enh_statistics/TS_enh_TF_targetgene_3


#找到enhancer与CpG位置overlap的并输出为xls表格：
awk 'BEGIN{FS="\t";OFS="\t";}{ if(NR>1){print $1,$2,$3,$4,$5,$6,$7} }' /home/ding/all_cmd/script/enh_statistics/enh_all_count_spe.bed|sortBed -i - |closest-features --closest  --delim '\t' --dist  - /media/ding/000B49000006264C/eRNA_project/CPG/CPG.bed|awk 'BEGIN{FS="\t";OFS="\t"}{ if($NF==0){print $0}  }' - >/home/ding/all_cmd/script/enh_statistics/CPG_enh_dist0

#然后使用deal_format.R生成xlsx表：


#做CpG岛在增强子上的富集程度:
cd /home/ding/all_cmd/script/
#1.获得CpG富集的bed->bigwig文件：

awk 'BEGIN{OFS="\t"}{ print $8,$9,$10,$4"_CGi" }' /home/ding/all_cmd/script/enh_statistics/CPG_enh_dist0|sortBed -i ->./CGI0.bed
awk 'BEGIN{OFS="\t"}{ print $1,$2,$3,$4 }' /home/ding/all_cmd/script/enh_statistics/CPG_enh_dist0|sortBed -i - >./enh_CGI0.bed

#只选择落在增强子内的CpG岛的部分作为bigwig文件：
bedops --intersect ./CGI0.bed ./enh_CGI0.bed |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"CGI_enh_"NR}' - >./CGI0_enh_intersect.bed

bedToBam -i ./CGI0_enh_intersect.bed -g /home/ding/all_cmd/hg19.chrom_24.sizes >/home/ding/all_cmd/script/enh_statistics/CPG_enh_dist0.bam

samtools index /home/ding/all_cmd/script/enh_statistics/CPG_enh_dist0.bam /home/ding/all_cmd/script/enh_statistics/CPG_enh_dist0.bam.bai

cd /home/ding/all_cmd/script
bamCoverage  --binSize 1 --numberOfProcessors max  -b /home/ding/all_cmd/script/enh_statistics/CPG_enh_dist0.bam -o ./CPG_enh_dist0.bw



#使用bwtool获得CGI在增强子的覆盖：
cd /home/ding/all_cmd/script
bwtool aggregate 2000:2000 ./enh_CGI0.bed ./CPG_enh_dist0.bw -fill=0 /dev/stdout >./enh_statistics/CPG_enh_dist0.mean_signal 
#不需要计算平均信号，所以每个位置都乘982


#删除中间文件：
rm ./CGI0.bed ./enh_CGI0.bed ./CGI0_enh_intersect.bed /home/ding/all_cmd/script/enh_statistics/CPG_enh_dist0.bam 

################################
#获得10个组织每个组织最特异的转录因子的调控网络：
#在TF_TS_score_heatmap.R中查看每个组织最特异的转录因子，结果发现特异性转录因子的阈值是小于-40.2，而Adrenal中的TSPV却大于阈值，也就是说不是组织特异的转录因子，结果如下：
gedit /home/ding/all_cmd/script/enh_statistics/TS_TF_each_tissue



#############
#TS_enh_TF_targetgene_3.2 里是所有转录因子调控TS enh的文件，接下来获得32个TS TF调控的 TS enh等关系：
#根据 TS_32TF_each_9tissue 文件的第一列：
cd ~/all_cmd/script/enh_statistics
awk 'BEGIN{OFS="\t"}{ if(NR>1&&NR==FNR){ a[$1]=$2 } if(NR>FNR && $3 in a){ print $0 } }' TS_32TF_each_9tissue TS_enh_TF_targetgene_3.2 >TS_enh_TF_targetgene_3.2_32






















