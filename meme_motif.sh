#先获得每个组织的enh的bed文件：
#只针对中心100bp的序列：
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if(NR>1&&$8==tissue){printf "%s\t%d\t%d\n",$1,($2+$3)/2-100,($2+$3)/2+100} }' enh_all_count_spe_melt|sortBed -i - >./meme_motif/$i"_enh.bed"
echo $i"完成1"
done

#非编码基因区域作为HOMER的 background:
cd /home/ding/all_cmd/script/enh_statistics
#randomBed -seed 2 -n 100000000 -l 400 -g ../../hg19.chrom_24.sizes |sortBed -i ->./meme_random.bed
#fastaFromBed -fi /home/ding/hg19/hg19_index.fa -bed ./meme_random.bed -fo ./meme_random.fasta
#获得基因外bed和#获得基因外bed的fasta：
awk 'BEGIN{OFS="\t"}{print $1,1,$2}'  ../../hg19.chrom_24.sizes|sortBed -i - |bedops --difference ../../hg19.chrom_24.bed /media/ding/000B49000006264C/eRNA_project/gencode/protein_gene_sorted.bed | fastaFromBed -fi /home/ding/hg19/hg19_index.fa -bed - -fo /media/ding/000B49000006264C/eRNA_project/HOMER/background.fasta


i=Adrenal
#针对Adrenal单独生成：
cd /home/ding/all_cmd/script/enh_statistics
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if(NR>1&&$8==tissue){printf "%s\t%d\t%d\n",$1,($2+$3)/2-100,($2+$3)/2+100} }' enh_all_count_spe_melt|sortBed -i - >./meme_motif/$i"_enh.bed"
cd /home/ding/all_cmd/script/enh_statistics/meme_motif
fastaFromBed -fi /home/ding/hg19/hg19_index.fa -bed $i"_enh.bed" -fo $i"_enh.fasta"


#从UCSC获得每个组织enh的fasta文件
#突然发现可以使用bedtools的fastaFromBed来提取bed为fasta:
cd /home/ding/all_cmd/script/enh_statistics/meme_motif
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  fastaFromBed -fi /home/ding/hg19/hg19_index.fa -bed $i"_enh.bed" -fo $i"_enh.fasta"
echo $i"完成2"
done



#本地化meme并获得motif:
cd /home/ding/all_cmd/script/enh_statistics/meme_motif
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  fasta_maxsize=`wc -c $i"_enh.fasta"|awk '{print $1}' -`
  echo $fasta_maxsize
  
  cd /home/ding/all_cmd/script/enh_statistics/meme_motif
  meme ./$i"_enh.fasta" -maxsize 10000000000000 -w 20 -maxw 20  -nmotifs 1 -evt 0.01  -mod oops -dna -oc ./"meme_"$i
  
  echo $i"完成3"
done

#meme ./Adrenal_enh.fasta   -maxsize 10000000 -mod oops -dna -oc ./meme_Adrenal3

./configure --prefix=$HOME/meme2 --with-url=http://meme-suite.org --enable-build-libxml2 --enable-build-libxslt --enable-serial


#Adrenal:
http://csbl.bmb.uga.edu/DMINDA/motif_annotation_prediction.php?jobid=2016111611226f
Yon can remember your jobid 2016111611226f

#这个竟然没有找到，只能重新划定中心范围。
http://csbl.bmb.uga.edu/DMINDA/motif_annotation_prediction.php?jobid=2016112073940f
Yon can remember your jobid 2016112073940f

http://csbl.bmb.uga.edu/DMINDA/motif_annotation_prediction.php?jobid=2016112060658f
Yon can remember your jobid 2016112060658f

#上下游300bp:
Link: http://csbl.bmb.uga.edu/DMINDA/motif_annotation_prediction.php?jobid=2016120603702f
Yon can remember your jobid 2016120603702f

#上下游50bp的：
Link: http://csbl.bmb.uga.edu/DMINDA/motif_annotation_prediction.php?jobid=2016120681136f
Yon can remember your jobid 2016120681136f

#Brain:
http://csbl.bmb.uga.edu/DMINDA/motif_annotation_prediction.php?jobid=2016111611638f
Yon can remember your jobid 2016111611638f
11 	0.00e+00 	1581 
#Breast 
http://csbl.bmb.uga.edu/DMINDA/motif_annotation_prediction.php?jobid=2016111612529f
Yon can remember your jobid 2016111612529f
10 	1.75e-229 	431 
#Heart
 http://csbl.bmb.uga.edu/DMINDA/motif_annotation_prediction.php?jobid=2016111613109f
Yon can remember your jobid 2016111613109f
12 	3.93e-248 	858 
#Liver
http://csbl.bmb.uga.edu/DMINDA/motif_annotation_prediction.php?jobid=2016111614042f
Yon can remember your jobid 2016111614042f
11    3.48e-229 	4071 
#Lung 
http://csbl.bmb.uga.edu/DMINDA/motif_annotation_prediction.php?jobid=2016111615048f
Yon can remember your jobid 2016111615048f
10 	5.69e-256 	1513
#Ovary 
http://csbl.bmb.uga.edu/DMINDA/motif_annotation_prediction.php?jobid=2016111615152f
Yon can remember your jobid 2016111615152f
11 	5.47e-73 	126 

#Placenta 
http://csbl.bmb.uga.edu/DMINDA/motif_annotation_prediction.php?jobid=2016111615546f
Yon can remember your jobid 2016111615546f
11 	0.00e+00 	994 

#SkeletalMuscle 
http://csbl.bmb.uga.edu/DMINDA/motif_annotation_prediction.php?jobid=2016111620317f
Yon can remember your jobid 2016111620317f
10 	7.26e-169 	548 
#Kidney
 http://csbl.bmb.uga.edu/DMINDA/motif_annotation_prediction.php?jobid=2016111620829f
Yon can remember your jobid 2016111620829f
12 	0.00e+00 	633 




#使用R包的seqLogo进行寻找每个组织的sequence logo:
#首先获得每个位置的每个碱基的出现频率。
cd /home/ding/all_cmd/script/enh_statistics/meme_motif
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  line_sum=`wc -l ./$i"_enh.fasta"|awk '{print $1}' -`
  awk  -v tissue=$i -v line_sum=$line_sum 'BEGIN{OFS="\t";FS=""}{ if(NR%2==0){  } }' ./$i"_enh.fasta" |head
  echo $i"完成3" 
done

#HOMER:
cd /media/ding/000B49000006264C/eRNA_project/HOMER
/home/ding/tool/HOMER/bin/findMotifs.pl  ../meme_motif/Adrenal_enh.fasta hg19 ./Adrenal -fasta /home/ding/hg19/hg19_index.fa

#HOMER:
cd /media/ding/000B49000006264C/eRNA_project/HOMER
/home/ding/tool/HOMER/bin/findMotifs.pl  ../meme_motif/Brain_enh.fasta hg19 ./Adrenal -fasta /home/ding/hg19/hg19_index.fa



#分两个进程：1
:>/media/ding/000B49000006264C/eRNA_project/HOMER/logo.txt
cd /home/ding/all_cmd/script/enh_statistics/meme_motif
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  #这个直接就可以：date +"%F %T"
  #记录起始时间$date_start并转换为时间$time_start：
  date_start=`date +"%Y-%m-%d %H:%M:%S"`
  time_start=`date +%s -d "$date_start"`
  
  #程序部分：
  cd /media/ding/000B49000006264C/eRNA_project/HOMER
  /home/ding/tool/HOMER/bin/findMotifs.pl  /home/ding/all_cmd/script/enh_statistics/meme_motif/$i"_enh_conservation.fasta"  fasta  ./$i -fastaBg  /home/ding/hg19/hg19_index.fa -len 18,20 -p 2 -noweight -chopify -bits
  echo $i"完成1" 
  
  #记录终止时间$date_end并转换为时间$time_end：
  date_end=`date +"%F %T"`
  time_end=`date +%s -d "$date_end"`
  
  #记录程序运行了多长时间：
  let time_range_hour=(time_end-time_start)/60/60
  let time_range_minute=(time_end-time_start)/60%60
  echo -e "组织:"$i"\t开始时间:"$date_start"\t结束时间:$date_end\t运行了"$time_range_hour"小时"$time_range_minute"分钟" >>./logo.txt
done


#####################################
#只能靠其他思路了：
#获得每个组织增强子保守性最高的前1000条序列：
#####################################

#先获得每个组织的enh的bed文件：
#只针对中心100bp的序列：
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if(NR>1&&$8==tissue){printf "%s\t%d\t%d\n",$1,($2+$3)/2-100,($2+$3)/2+100} }' enh_all_count_spe_melt|sortBed -i - >./meme_motif/$i"_enh.bed"
echo $i"完成1"
done

#获得保守性信号，sum(mean)最高的前2000条，将NA替换为0
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  bwtool summary ./meme_motif/$i"_enh.bed" /media/ding/000B49000006264C/eRNA_project/phastCons46/phastCons46.bw -with-sum /dev/stdout -fill=0 |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$NF}' -|sort  -k 4 -n -r - |head -n 2000|sortBed -i - >./meme_motif/$i"_enh_conservation.bed"
echo $i"完成2"
done

#使用bedtools的fastaFromBed来提取bed为fasta:
cd /home/ding/all_cmd/script/enh_statistics/meme_motif
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  fastaFromBed -fi /home/ding/hg19/hg19_index.fa -bed $i"_enh_conservation.bed" -fo $i"_enh_conservation.fasta"
echo $i"完成3"
done
###########################################



###########################################
#获得DHS信号最高的前2000条进行寻找motif:

#先获得每个组织的enh的bed文件：
#只针对中心100bp的序列：
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if(NR>1&&$8==tissue){printf "%s\t%d\t%d\n",$1,($2+$3)/2-100,($2+$3)/2+100} }' enh_all_count_spe_melt|sortBed -i - >./meme_motif/$i"_enh.bed"
echo $i"完成1"
done


#获得DHS信号，sum(mean)最高的前3000条，将NA替换为0 -fill=0
#获得DHS信号最高的前50%条：
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  bwtool summary ./meme_motif/$i"_enh.bed" /media/ding/000B49000006264C/eRNA_project/DNase/$i/new_DHS_TPM.bw -with-sum /dev/stdout -fill=0 |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$NF}' -|sort -k 4 -n -r  |head -n 2000|sortBed -i - >./meme_motif/$i"_enh_DHS.bed"
 
  #获得每个组织增强子的数目的50%： 
#  enh_num=`wc -l ./meme_motif/$i"_enh.bed"|awk '{print $1}' -`
#  let  enh_num_meme=$enh_num/2
#  bwtool summary ./meme_motif/$i"_enh.bed" /media/ding/000B49000006264C/eRNA_project/DNase/$i/new_DHS_TPM.bw -with-sum /dev/stdout -fill=0 |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$NF}' -|sort  -k 4 -n -r - |head -n $enh_num_meme|sortBed -i ->./meme_motif/$i"_enh_DHS.bed"
  
echo $i"完成2"
done


#使用bedtools的fastaFromBed来提取bed为fasta:
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  fastaFromBed -fi /home/ding/hg19/hg19_index.fa -bed ./meme_motif/$i"_enh_DHS.bed" -fo ./meme_motif/$i"_enh_DHS.fasta"
echo $i"完成3"
done

###########################################

#参数变量：
meme_args=" -revcomp -maxsize 100000000000 -nmotifs 1 -minw 12 -evt 0.01  -mod oops -dna"

#设置函数：
function do_meme(){
  if [ ${1} != '' ] 
  then
    tissue=${1}
    echo $tissue
    cd /home/ding/all_cmd/script/enh_statistics/meme_motif 
    echo $meme_args|xargs meme -oc ./meme_${tissue} ./${tissue}_enh_DHS.fasta
  fi
}


#开始执行：
do_meme Adrenal 
do_meme Brain
do_meme Breast
do_meme Heart
do_meme Liver
do_meme Lung
do_meme Ovary
do_meme Placenta
do_meme SkeletalMuscle
do_meme Kidney



















