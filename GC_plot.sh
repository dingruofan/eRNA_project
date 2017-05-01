#!/bin/bash
#使用的是没有strand的TSS：
cd /home/ding/all_cmd/script/
n=`wc -l ./enh_statistics/enh_all_count|awk '{print $1-1}' - `
randomBed -l 200 -n $n -seed 2 -g ../hg19.chrom_24.sizes|sortBed -i - |bedops --not-element-of 1 - /media/ding/000B49000006264C/eRNA_project/gencode/protein_coding_gene.bed |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"random_"NR}' ->./enh_statistics/random_all_enh.bed

bwtool agg 2000:2000 ./enh_statistics/random_all_enh.bed,./enh_statistics/enh_bid_cluster.bed,./enh_statistics/enh_unbid_cluster.bed,/media/ding/000B49000006264C/eRNA_project/gencode/protein_gene2.ss,./enh_statistics/no_eRNA_enh.bed /media/ding/000B49000006264C/eRNA_project/GC/hg19.gc5Base.bw  /dev/stdout > ./enh_statistics/GC_result.mean_signal











#############################

#每个组织中的每个增强子序列的GC含量：
########

tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  #先获得每个组织的enh中心100bp的bed文件：
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if(NR>1&&$8==tissue){printf "%s\t%d\t%d\t%s\n",$1,($2+$3)/2-100,($2+$3)/2+100,$4} }' ./enh_all_count_spe_melt|sortBed -i - >./GC/bed/$i"_enh.bed"
  echo $i"完成1"

  # 然后获得每个bed文件的fasta文件：
  fastaFromBed -fi /home/ding/hg19/hg19_index.fa -bed ./GC/bed/$i"_enh.bed" -fo ./GC/fasta/$i"_enh.fasta"
  echo $i"完成2"

  #然后使用tool的脚本获得每个enh所对应序列的GC含量：
   get_gc_content.pl ./GC/fasta/$i"_enh.fasta" #获得一个gc_out.txt结果文件：
   #使用awk进行处理并删除结果文件：
   awk -v tissue="$i" 'BEGIN{OFS="\t";FS="\t";print "ID","GC_percent","tissue"}{ if(NR>1){print $1,$2,tissue} }' ./gc_out.txt >./GC/GC_percent/$i"_enh_GC" 
   rm ./gc_out.txt
   echo $i"完成3"

   #添加增强子的名称到结果文件的最后一列：
   awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[">"$1":"$2"-"$3]=$4; } if(NR>FNR&&FNR==1){print $0,"enh_name"} if(NR>FNR&&FNR>1){print $0,a[$1]} }' ./GC/bed/$i"_enh.bed" ./GC/GC_percent/$i"_enh_GC" >./GC/GC_percent/a
   awk '{print $0}' ./GC/GC_percent/a >./GC/GC_percent/$i"_enh_GC"
   rm ./GC/GC_percent/a
done

cd /home/ding/all_cmd/script/enh_statistics
#使用上面生成的：./enh_statistics/random_all_enh.bed
#获得20000个基因外随机位点：
randomBed -l 400 -n 90000 -seed 2 -g /home/ding/all_cmd/hg19.chrom_24.sizes|sortBed -i - |bedops --not-element-of 1 - /media/ding/000B49000006264C/eRNA_project/gencode/protein_coding_gene.bed |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"random_"NR}' - >./GC/bed/random.bed
echo "随机位点bed完成~"

# 然后获得每个bed文件的fasta文件：
fastaFromBed -fi /home/ding/hg19/hg19_index.fa -bed ./GC/bed/random.bed -fo ./GC/fasta/random.fasta
echo "随机位点fasta完成~"

#然后使用tool的脚本获得每个enh所对应序列的GC含量：
get_gc_content.pl ./GC/fasta/random.fasta #获得一个gc_out.txt结果文件：
#使用awk进行处理并删除结果文件：
awk -v tissue="random" 'BEGIN{OFS="\t";FS="\t";print "ID","GC_percent","tissue"}{ if(NR>1){print $1,$2,tissue} }' ./gc_out.txt >./GC/GC_percent/random_GC
rm ./gc_out.txt
#添加增强子的名称到结果文件的最后一列：
awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[">"$1":"$2"-"$3]=$4; } if(NR>FNR&&FNR==1){print $0,"enh_name"} if(NR>FNR&&FNR>1){print $0,a[$1]} }' ./GC/bed/random.bed ./GC/GC_percent/random_GC>./GC/GC_percent/a
awk '{print $0}' ./GC/GC_percent/a >./GC/GC_percent/random_GC
rm ./GC/GC_percent/a
echo "随机位点GC统计完成3"     


#将十个组织以及随机位点的GC合并：
cd /home/ding/all_cmd/script/enh_statistics
#添加列名并创建文件：
awk '{print $0}' ./GC/GC_percent/random_GC >./GC/GC_all.tsv
#开始添加内容：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  awk '{if(NR>1){print $0}}' ./GC/GC_percent/$i"_enh_GC" >>./GC/GC_all.tsv
  echo $i"完成4"
done


########



##################
###统计每个增强子的GC并求平均值：
cd /home/ding/all_cmd/script/enh_statistics
fastaFromBed -fi /home/ding/hg19/hg19_index.fa -bed ./enh_all_cluster.bed -fo ./GC/fasta/all_enh.fasta
get_gc_content.pl ./GC/fasta/all_enh.fasta #获得一个gc_out.txt结果文件：
#使用awk进行处理并删除结果文件：
awk -v tissue="all_enh" 'BEGIN{OFS="\t";FS="\t";print "ID","GC_percent","all_enh"}{ if(NR>1){print $1,$2,tissue} }' ./gc_out.txt >./GC/GC_percent/all_enh_GC
rm ./gc_out.txt
#获得平均GC含量：
awk 'BEGIN{sum_GC=0;nrow=0}{if(NR>1){sum_GC=sum_GC+$2;nrow=NR}}END{print sum_GC,nrow,sum_GC/nrow}' ./GC/GC_percent/all_enh_GC
#结果：总GC：2.53658e+06 增强子数目：53925  平均每个增强子GC含量：47.039
echo "all_enh GC统计完成3"     

##########################















