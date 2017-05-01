#!/bin/bash

#######################################################################
#根input归一化后的组蛋白修饰信号计算组织特异性：
#先使用bwtool提取enh中心附近2000bp内的组蛋白信号平均值：
# 
# 
# #先提取enh中心附近2000bp内区域：
# cd /home/ding/all_cmd/script/enh_statistics
# awk 'BEGIN{FS="\t";OFS="\t"}{ printf "%s\t%.f\t%.f\t%s\n",$1,($2+$3)/2-2000,($2+$3)/2+2000,$4 }' ./enh_all_cluster.bed >./enh_all_cluster_200.bed
# 
# #先创建一个平均信号值矩阵：
# cd /home/ding/all_cmd/script/enh_statistics
# awk 'BEGIN{FS="\t";OFS="\t";print "chr","chr_start","chr_end","enh"}{ print $0 }' ./enh_all_cluster.bed >./enh_all_cluster_histone_sig.matrix
# #计算并填充矩阵：
# tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
# #tissue=(Adrenal)
# for i in ${tissue[@]}
# do 
# cd /media/ding/000B49000006264C/eRNA_project/histone/$i
# #第八列为平均值：
# bwtool summary -fill=0 /home/ding/all_cmd/script/enh_statistics/enh_all_cluster_200.bed ./bg/H3K4me1.bw /dev/stdout  |awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$1"\t"$2"\t"$3]=$8} if(NR>FNR){print $4,a[$1"\t"$2"\t"$3]} }' - /home/ding/all_cmd/script/enh_statistics/enh_all_cluster_200.bed |awk -v tissue=$i 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$1]=$2} if(NR>FNR&&FNR==1){print $0,tissue} if(NR>FNR&&FNR>1){print $0,a[$4]} }' - /home/ding/all_cmd/script/enh_statistics/enh_all_cluster_histone_sig.matrix >/home/ding/all_cmd/script/enh_statistics/a
# cp /home/ding/all_cmd/script/enh_statistics/a  /home/ding/all_cmd/script/enh_statistics/enh_all_cluster_histone_sig.matrix
# rm /home/ding/all_cmd/script/enh_statistics/a 
# echo $i"完成"
# done 

#根据公式计算组织特异性值，在R下计算更方便：
cd /home/ding/all_cmd/script
Rscript ./enh_ts.R

#####################################
#生成了enh特异性注释文件。./enh_statistics/enh_spe
#现有特异性的TF文件：./enh_statistics/TF_spe

#寻找特异性的增强子与特异性的TF之间的对应关系：

#先找出这些转录因子的结合位点数据，保留TFBS，TF_name,TF_score的bed格式：
cd /home/ding/all_cmd/script
awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR>1&&NR==FNR){a[$1]=$NF} if(NR>FNR&&$4 in a){print $1,$2,$3,$4,a[$4]} }' ./enh_statistics/TF_spe /media/ding/000B49000006264C/eRNA_project/TFBS/TFBS_sorted.bed |sortBed -i - >./enh_statistics/TF_spe_TFBS

#每个enh和组织对应起来。
awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==1){for(i=1;i<=NF;i++){a[i]=$i}} if(NR>1){printf "%s\t%s\t%s\t%s\t",$1,$2,$3,$4;isFirst="TURE";for(j=5;j<=NF-2;j++){if($j==1){ if(isFirst=="FALSE"){printf ",%s",a[j]}  if(isFirst=="TURE"){printf "%s",a[j];isFirst="FALSE"}  }}} print "" }' ./enh_statistics/enh_all_count| sortBed -i -  >./enh_statistics/enh_all_count_tissue.bed

#特异性enh附近2000bp内特异性的TF：
cd /home/ding/all_cmd/script
awk '(NR>1&&$7=="spe"){print}' ./enh_statistics/enh_spe |sortBed -i - | closest-features --closest  --delim '\t' --dist ./enh_statistics/TF_spe_TFBS -  |awk '{if($NF>=-2000&&$NF<=2000){print }}' - |awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$4]=$5} if(NR>FNR){print $0,a[$9]} }' ./enh_statistics/enh_all_count_tissue.bed - >./enh_statistics/a
#整理格式：一个enh一行保留多个逗号分割的TF,用R整理生成文件：./enh_statistics/spe_TF_enh_tissue_1
cd /home/ding/all_cmd/script
Rscript ./enh_ts2.R
rm ./enh_statistics/a

#再用awk整理格式：
cd /home/ding/all_cmd/script
#整理为TF和tissue都分隔到单行：
awk 'BEGIN{FS="\t";OFS="\t"}{ split($4,a,","); for(i=1;i<=length(a);i++){print $1,$2,$3,a[i]} }' ./enh_statistics/spe_TF_enh_tissue_1 |awk 'BEGIN{FS="\t";OFS="\t"}{ split($3,a,","); for(i=1;i<=length(a);i++){print $1,$2,$4,a[i]} }' - >./enh_statistics/spe_TF_enh_tissue_2

#整理tissue到单行，TF还在一行：
awk 'BEGIN{FS="\t";OFS="\t"}{ split($3,a,","); for(i=1;i<=length(a);i++){print $1,$2,$4,a[i]} }' ./enh_statistics/spe_TF_enh_tissue_1 >./enh_statistics/spe_TF_enh_tissue_3

gedit ./enh_statistics/spe_TF_enh_tissue_1
gedit ./enh_statistics/spe_TF_enh_tissue_2


#对所有的增强子进行全面注释：
#得先运行enh_ts.R，得到enh_all_count_spe:
cd /home/ding/all_cmd/script/enh_statistics
awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==1){for(i=1;i<=NF;i++){a[i]=$i}} if(NR>1){printf "%s\t%s\t%s\t%s\t%s\t%s\t",$1,$2,$3,$4,$16,$17;isFirst="TURE";for(j=5;j<=14;j++){if($j==1){ if(isFirst=="FALSE"){printf ",%s",a[j]}  if(isFirst=="TURE"){printf "%s",a[j];isFirst="FALSE"}  }}} print "" }' ./enh_all_count_spe| sortBed -i -  >./enh_all_count_spe.bed
#然后使用/home/ding/all_cmd/script/deal_format.R输出为xlsx

#获得每个组织的增强子及其TFBS/TF exp
#处理：./enh_find/enh_bid_TF.bed ./enh_find/enh_no_eRNA_named_TFBS_exp ./enh_find_1/enh_unbid_TF.bed
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  #先整理格式统一列，然后保存到一个文件
  awk 'BEGIN{OFS="\t"}{ print $1,$2,$3,$4,"2D-eRNA",$6,$7,$8,$9,$10,$11 }' ./enh_find/enh_bid_TF.bed >./enh_find/a1
  awk 'BEGIN{OFS="\t"}{ print $1,$2,$3,$4,"1D-eRNA",$6,$7,$8,$9,$10,$11 }' ./enh_find_1/enh_unbid_TF.bed >./enh_find/a2
  awk 'BEGIN{OFS="\t"}{ print $1,$2,$3,$4,"no eRNAs",".",$5,$6,$7,$8,$10 }' ./enh_find/enh_no_eRNA_named_TFBS_exp >./enh_find/a3
  
  awk 'BEGIN{OFS="\t"}{ print }' ./enh_find/a1 ./enh_find/a2 |awk 'BEGIN{OFS="\t"}{ print }' - ./enh_find/a3 >./enh_find/enh_all_TF.bed
  
  rm ./enh_find/a1 ./enh_find/a2 ./enh_find/a3
  echo $i"_完成"
done

#然后使用/home/ding/all_cmd/script/deal_format.R输出为xlsx









