#!/bin/bash

#######################################################################
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
tissue=(Adrenal)
for i in ${tissue[@]}
do 
  
  cd /media/ding/000B49000006264C/eRNA_project/fantom/$i"_CAGE"
  #得到正链ctss的bam文件：
  awk 'BEGIN{FS="\t";OFS="\t"}{ if($NF=="+"){for(i=1;i<=$5;i++){print $1,$2,$3,"ctss_"NR,1,$NF}} }' ctss.bed  |bedToBam -i - -g /home/ding/all_cmd/hg19.chrom_24.sizes >./ctss_pos.bam
  samtools sort ./ctss_pos.bam -o ./ctss_pos_sorted.bam
  samtools index ./ctss_pos_sorted.bam ./ctss_pos_sorted.bam.bai
  #得到负链ctss的bam文件：
  awk 'BEGIN{FS="\t";OFS="\t"}{ if($NF=="-"){for(i=1;i<=$5;i++){print $1,$2,$3,"ctss_"NR,1,$NF}} }' ctss.bed  |bedToBam -i - -g /home/ding/all_cmd/hg19.chrom_24.sizes >./ctss_neg.bam
  samtools sort ./ctss_neg.bam -o ./ctss_neg_sorted.bam
  samtools index ./ctss_neg_sorted.bam ./ctss_neg_sorted.bam.bai
  #信号和用于归一化：
  #正链的ctss信号和：sum_sig_pos=`awk 'BEGIN{}{ if($NF=="+"){for(i=1;i<=$5;i++){print }} }' ctss.bed |wc -l`
  #负链的ctss信号和：sum_sig_neg=`awk 'BEGIN{}{ if($NF=="-"){for(i=1;i<=$5;i++){print }} }' ctss.bed |wc -l`
  
  #使用deeptools的bamcoverage 可以直接将bam文件转为信号文件bw/bg，而且是直接归一化的，妈的，histone那里走弯路了！草！，直接bamCoverage就完事了。。甚至直接--filterRNAstrand就不用上面的正负链分类了！
  bamCoverage --scaleFactor 1.0  --binSize 1 --numberOfProcessors max  -b ./ctss_pos_sorted.bam -o ./pos.bw 
  bamCoverage --scaleFactor 1.0  --binSize 1 --numberOfProcessors max  -b ./ctss_neg_sorted.bam -o ./neg.bw 
  echo $i"完成1"
done

#Liver
#bamCoverage --scaleFactor 1.0 --filterRNAstrand forward --binSize 100  --numberOfProcessors max  -b ./liver%2c%20adult%2c%20pool1.CNhs10624.10018-101C9.hg19.nobarcode.bam -o ./pos.bw 
#bamCoverage --scaleFactor 1.0 --filterRNAstrand reverse --binSize 100  --numberOfProcessors max  -b ./liver%2c%20adult%2c%20pool1.CNhs10624.10018-101C9.hg19.nobarcode.bam -o ./neg.bw 


#生成每个组织pos和neg的ctss文件：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
tissue=(Adrenal)
for i in ${tissue[@]}
do 
  cd /media/ding/000B49000006264C/eRNA_project/fantom/$i"_CAGE"
  awk 'BEGIN{FS="\t";OFS="\t"}{ if($NF=="+"){for(i=1;i<=$5;i++){print $1,$2,$3,"ctss_"NR,1,$NF}} }' ctss.bed |sortBed -i - >ctss_pos.bed
  awk 'BEGIN{FS="\t";OFS="\t"}{ if($NF=="-"){for(i=1;i<=$5;i++){print $1,$2,$3,"ctss_"NR,1,$NF}} }' ctss.bed |sortBed -i - >ctss_neg.bed
  echo $i"完成2"
done



#每个组织的bid/unbid合并后，用bwtool统计其-2000:2000的ctss分布并做图：

tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  #Breast Liver Lung
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  #合并bid/unbid的enh到临时文件a.bed
  awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$4}' ./enh_find/enh_bid.bed ./enh_find_1/enh_unbid.bed  |awk 'BEGIN{OFS="\t"}{printf "%s\t%d\t%d\t%s\n",$1,($2+$3)/2-0,($2+$3)/2+1,$4}' - |sortBed -i - >./a.bed
  bwtool agg 5000:5000 ./a.bed ../../fantom/$i"_CAGE"/pos.bw,../../fantom/$i"_CAGE"/neg.bw  /dev/stdout >./result.mean_signal

  #用R做图
  #now_path=`pwd`
  #Rscript /home/ding/all_cmd/script/enh_ctss_enrich.R $now_path
  #feh ./result.mean_signal.png
  #rm ./a.bed ./result.mean_signal ./result.mean_signal.png
  echo $i"完成3"
done

#重新测试Liver的：
# 

  echo $i
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  # ./enh_find/enh_bid.bed ./enh_find_1/enh_unbid.bed ./enh_find/enh_no_eRNA_named.bed
  awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$4}' ./enh_find/result.bed |awk 'BEGIN{OFS="\t"}{printf "%s\t%d\t%d\n",$1,($2+$3)/2-0,($2+$3)/2+1}' - |sortBed -i ->./a.bed
  closest-features --closest --dist --delim "\t"   ../../fantom/$i"_CAGE"/ctss_pos.bed ./a.bed |awk '{if($NF>=-500&&$NF<=500)print $NF}' - |sort - |uniq -c - |awk '{print $2"\t"$1}' - >./pos.dist
  closest-features --closest --dist --delim "\t"   ../../fantom/$i"_CAGE"/ctss_neg.bed ./a.bed |awk '{if($NF>=-500&&$NF<=500)print $NF}' - |sort - |uniq -c - |awk '{print $2"\t"$1}' - >./neg.dist
# 





###########################################################
#合并所有组织的正链和负链的ctss，并转为bw文件：
:> /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/ctss_pos.bed
:> /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/ctss_neg.bed
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  cd /media/ding/000B49000006264C/eRNA_project/fantom/
  cat $i"_CAGE"/ctss_pos.bed>>./All_tissue_CAGE/ctss_pos.bed  
  cat $i"_CAGE"/ctss_neg.bed>>./All_tissue_CAGE/ctss_neg.bed  
  echo $i"完成1"
done

#将合并后的bed转为bam后再转为bw:
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/
#转为bam:
  bedToBam -i ./ctss_pos.bed  -g /home/ding/all_cmd/hg19.chrom_24.sizes >./ctss_pos.bam
  bedToBam -i ./ctss_neg.bed  -g /home/ding/all_cmd/hg19.chrom_24.sizes >./ctss_neg.bam
#对bam排序并建立索引：
  samtools sort ./ctss_pos.bam -o ./ctss_pos_sorted.bam
  samtools index ./ctss_pos_sorted.bam ./ctss_pos_sorted.bam.bai
  samtools sort ./ctss_neg.bam -o ./ctss_neg_sorted.bam
  samtools index ./ctss_neg_sorted.bam ./ctss_neg_sorted.bam.bai
#转为bw文件：
 bamCoverage --scaleFactor 1.0  --binSize 1 --numberOfProcessors max  -b ./ctss_pos_sorted.bam -o ./pos.bw 
 bamCoverage --scaleFactor 1.0  --binSize 1 --numberOfProcessors max  -b ./ctss_neg_sorted.bam -o ./neg.bw 
#附近信号：
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/
bwtool agg 5000:5000 /home/ding/all_cmd/script/enh_statistics/enh_all_cluster.bed ./pos.bw,./neg.bw  /dev/stdout >./result.mean_signal


#获得所有增强子与所有gene的位置关系：
cd /home/ding/all_cmd/script/enh_statistics
closest-features --closest --dist --delim "\t" /media/ding/000B49000006264C/eRNA_project/gencode/all_gene_sorted.bed  ./enh_all_count_spe.bed >./enh_allgene_dist






