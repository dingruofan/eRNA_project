#!/bin/bash

#每个组织增强子中心100bp内GC含量的均值与其靶基因的表达量是否有相关性。




#首先获得每个增强子中心100bp内GC含量的均值
cd /home/ding/all_cmd/script/enh_statistics
awk 'BEGIN{OFS="\t"}{ if(NR>0){printf "%s\t%d\t%d\t%s\n",$1,($2+$3)/2-100,($2+$3)/2+100,$4} }' ./enh_all_cluster.bed |sortBed -i>./a
bwtool summary ./a /media/ding/000B49000006264C/eRNA_project/GC/hg19.gc5Base.bw  -with-sum /dev/stdout | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$8}'>./b
#获得每个增强子中心100bp内的GC含量的均值。
awk 'BEGIN{OFS="\t";FS="\t"}{ if(NR==FNR){a[$1"\t"$2"\t"$3]=$4} if(NR>FNR){print a[$1"\t"$2"\t"$3],$4} }' ./a ./b >./c

#合并到文件enh_all_count_spe_melt的最后一列：
awk 'BEGIN{OFS="\t";FS="\t"}{ if(NR==FNR){a[$1]=$2} if(NR>FNR&&FNR==1){print $0,"GC_mean"} if(NR>FNR&&FNR>1){print $0,a[$4]} }' ./c ./enh_all_count_spe_melt >./enh_all_count_spe_melt_GC
rm ./a ./b ./c 




#从enh_exp_targetgene_exp获得各个组织的增强子的靶基因的表达量，
#注意:placenta的表达量为-1.
cd /home/ding/all_cmd/script/enh_statistics
#创建空文件用来反复>>。
:>./enh_GC_target_eRNA_exp
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 

  
  echo
done
#增强子有靶基因的输出，补充基因type和GC_mean。
  awk -v tissue=$i 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR&&NR>1){a[$4"\t"$8]=$7"\t"$9;}  if(NR>FNR&&($4"\t"$13 in a)){print $0,a[$4"\t"$13]} if(NR>FNR&&!($4"\t"$13 in a)){print $0,"random_type","random_GC"} }'  ./enh_all_count_spe_melt_GC ./enh_exp_targetgene_exp >./enh_GC_target_eRNA_exp



















