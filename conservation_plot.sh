#!/bin/bash
cd /home/ding/all_cmd/script/
n=`wc -l ./enh_statistics/enh_all_count|awk '{print $1-1}' - `
randomBed -l 200 -n $n -seed 2 -g ../hg19.chrom_24.sizes|sortBed -i - |bedops --not-element-of 1 - /media/ding/000B49000006264C/eRNA_project/gencode/protein_coding_gene.bed >./enh_statistics/random_all_enh.bed

bwtool agg 2000:2000 ./enh_statistics/random_all_enh.bed,./enh_statistics/enh_bid_cluster.bed,./enh_statistics/enh_unbid_cluster.bed,/media/ding/000B49000006264C/eRNA_project/gencode/protein_gene2.ss,./enh_statistics/no_eRNA_enh.bed /media/ding/000B49000006264C/eRNA_project/phastCons46/phastCons46.bw  /dev/stdout > ./enh_statistics/conservation_result.mean_signal




#############################

#每个组织中的每个增强子序列的conservation分值：
########

tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  #先获得每个组织的enh中心100bp的bed文件：
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if(NR>1&&$8==tissue){printf "%s\t%d\t%d\t%s\n",$1,($2+$3)/2-100,($2+$3)/2+100,$4} }' ./enh_all_count_spe_melt|sortBed -i - >./conservation/bed/$i"_enh.bed"
  echo $i"完成1"
  
  #通过bwtools获得每个增强子上的平均保守分数：
  bwtool summary ./conservation/bed/$i"_enh.bed" /media/ding/000B49000006264C/eRNA_project/phastCons46/phastCons46.bw -with-sum  -fill=0 /dev/stdout |awk -v tissue="$i" 'BEGIN{OFS="\t";FS="\t";print "ID","phastCons","tissue"}{print $1":"$2"-"$3,$8,tissue}' - >./conservation/conservation_percent/$i"_enh_conservation" 
   echo $i"完成2"
done


cd /home/ding/all_cmd/script/enh_statistics
#使用上面生成的：./enh_statistics/random_enh.bed
#获得20000个基因外随机位点：
randomBed -l 400 -n 90000 -seed 2 -g /home/ding/all_cmd/hg19.chrom_24.sizes|sortBed -i - |bedops --not-element-of 1 - /media/ding/000B49000006264C/eRNA_project/gencode/protein_coding_gene.bed |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"random_"NR}' - >./conservation/bed/random.bed
echo "随机位点bed完成~"

 bwtool summary ./conservation/bed/random.bed /media/ding/000B49000006264C/eRNA_project/phastCons46/phastCons46.bw -with-sum -fill=0 /dev/stdout  |awk 'BEGIN{OFS="\t";FS="\t";print "ID","phastCons","tissue"}{print $1":"$2"-"$3,$8,"random"}' - >./conservation/conservation_percent/random_conservation
echo "随机位点conservation统计完成2"     


#将十个组织以及随机位点的conservation合并：
cd /home/ding/all_cmd/script/enh_statistics
#添加列名并创建文件：
awk '{print $0}' ./conservation/conservation_percent/random_conservation >./conservation/conservation_all.tsv
#开始添加内容：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  awk '{if(NR>1){print $0}}' ./conservation/conservation_percent/$i"_enh_conservation" >>./conservation/conservation_all.tsv
  echo $i"完成4"
done
#添加列名：


########
#获得2D enh和1D enh的每个增强子的保守性
#先获得每个组织的enh中心100bp的bed文件：
cd /home/ding/all_cmd/script/enh_statistics
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$16=="bid"){printf "%s\t%d\t%d\t%s\n",$1,($2+$3)/2-100,($2+$3)/2+100,$4} }' ./enh_all_count_spe|sortBed -i - >./conservation/bed/bid_enh.bed
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$16=="unbid"){printf "%s\t%d\t%d\t%s\n",$1,($2+$3)/2-100,($2+$3)/2+100,$4} }' ./enh_all_count_spe|sortBed -i - >./conservation/bed/unbid_enh.bed
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$16=="no_eRNA"){printf "%s\t%d\t%d\t%s\n",$1,($2+$3)/2-100,($2+$3)/2+100,$4} }' ./enh_all_count_spe|sortBed -i - >./conservation/bed/no_eRNA_enh.bed
  
#通过bwtools获得每个增强子上的平均保守分数：
bwtool summary ./conservation/bed/bid_enh.bed  /media/ding/000B49000006264C/eRNA_project/phastCons46/phastCons46.bw -with-sum -fill=0 /dev/stdout |awk 'BEGIN{OFS="\t";FS="\t";print "ID","phastCons","type"}{print $1":"$2"-"$3,$8,"bid"}' - >./conservation/conservation_percent/bid_enh_conservation

bwtool summary ./conservation/bed/unbid_enh.bed  /media/ding/000B49000006264C/eRNA_project/phastCons46/phastCons46.bw -with-sum -fill=0 /dev/stdout |awk 'BEGIN{OFS="\t";FS="\t";print "ID","phastCons","type"}{print $1":"$2"-"$3,$8,"unbid"}' - >./conservation/conservation_percent/unbid_enh_conservation

bwtool summary ./conservation/bed/no_eRNA_enh.bed  /media/ding/000B49000006264C/eRNA_project/phastCons46/phastCons46.bw -with-sum -fill=0 /dev/stdout |awk 'BEGIN{OFS="\t";FS="\t";print "ID","phastCons","type"}{print $1":"$2"-"$3,$8,"no_eRNA"}' - >./conservation/conservation_percent/no_eRNA_enh_conservation


#然后使用R计算是否有显著差异：



########
#获得TS enh和UE enh的每个增强子的保守性
#先获得每个组织的enh中心100bp的bed文件：
cd /home/ding/all_cmd/script/enh_statistics
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$17=="spe"){printf "%s\t%d\t%d\t%s\n",$1,($2+$3)/2-100,($2+$3)/2+100,$4} }' ./enh_all_count_spe|sortBed -i - >./conservation/bed/TS_enh.bed
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$17=="other"){printf "%s\t%d\t%d\t%s\n",$1,($2+$3)/2-100,($2+$3)/2+100,$4} }' ./enh_all_count_spe|sortBed -i - >./conservation/bed/other_enh.bed
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$17=="uni"){printf "%s\t%d\t%d\t%s\n",$1,($2+$3)/2-100,($2+$3)/2+100,$4} }' ./enh_all_count_spe|sortBed -i - >./conservation/bed/UE_enh.bed
  
#通过bwtools获得每个增强子上的平均保守分数：
bwtool summary ./conservation/bed/TS_enh.bed  /media/ding/000B49000006264C/eRNA_project/phastCons46/phastCons46.bw -with-sum -fill=0 /dev/stdout |awk 'BEGIN{OFS="\t";FS="\t";print "ID","phastCons","type"}{print $1":"$2"-"$3,$8,"TS"}' - >./conservation/conservation_percent/TS_enh_conservation

bwtool summary ./conservation/bed/other_enh.bed  /media/ding/000B49000006264C/eRNA_project/phastCons46/phastCons46.bw -with-sum -fill=0 /dev/stdout |awk 'BEGIN{OFS="\t";FS="\t";print "ID","phastCons","type"}{print $1":"$2"-"$3,$8,"other"}' - >./conservation/conservation_percent/other_enh_conservation

bwtool summary ./conservation/bed/UE_enh.bed  /media/ding/000B49000006264C/eRNA_project/phastCons46/phastCons46.bw -with-sum -fill=0 /dev/stdout |awk 'BEGIN{OFS="\t";FS="\t";print "ID","phastCons","type"}{print $1":"$2"-"$3,$8,"uni"}' - >./conservation/conservation_percent/UE_enh_conservation


#然后使用R计算是否有显著差异：

#查看phastcons的保守性分值的范围：
cd /home/ding/all_cmd
sortBed -i ./hg19.chrom_24.bed >hg19.chrom_24_sorted.bed
bwtool summary ./hg19.chrom_24_sorted.bed /media/ding/000B49000006264C/eRNA_project/phastCons46/phastCons46.bw -header -fill=0 /dev/stdout >./script/enh_statistics/conservation/hg19.chrom_24_max_phastcons




