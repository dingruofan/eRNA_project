#!/bin/bash

#######################################################################

#获得每个组织的isspe_enh到文件./enh_find/enh_loc_type_isspe
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney )
#tissue=(Adrenal)
for i in ${tissue[@]}
do 
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  #首先获得每个组织的enh,通过处理 /home/ding/all_cmd/script/enh_statistics/enh_all_count_spe文件：
  awk -v tissue=$i 'BEGIN{OFS="\t";tissue_id=-1}{ if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if(NR>1&&$tissue_id==1){print $1,$2,$3,$4,$16,$17}  }' /home/ding/all_cmd/script/enh_statistics/enh_all_count_spe|sortBed -i - >./enh_find/enh_loc_type_isspe

  echo $i"完成1"
done

#使用enh_GO2.R分析：
#获得每个组织的enh的bed4格式到bed文件：
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
for i in ${tissue[@]}
do 
  awk -v tissue=$i 'BEGIN{OFS="\t"}{if($NF==tissue){print $1,$2,$3,$4}}' enh_all_count_spe_melt|sortBed -i - >./GO/bed/$i".bed4"
  echo $i"完成"
done


###获得每个组织的enh的靶基因，用来做GO分析：
#记得要去除掉random的：
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  awk -v tissue="$i" 'BEGIN{OFS="\t"}{ if($13==tissue&&$12!="random"){print $10} }' /home/ding/all_cmd/script/enh_statistics/enh_exp_targetgene_exp |sort -u >./target_gene/all_gene/$i"_target"
  echo $i" 完成target gene"
done

#合并十个组织的靶基因到一个文件：
cd /home/ding/all_cmd/script/enh_statistics
awk 'BEGIN{OFS="\t";print "gene_name","tissue"}' >./target_gene/all_gene/all_target
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  awk -v tissue="$i" '{print $0,tissue}' ./target_gene/all_gene/$i"_target" >>./target_gene/all_gene/all_target
done


#获得每个组织TS enh的靶基因：
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$4]=$6} if(NR>FNR){print $0,a[$4]} }' ./enh_all_count_spe.bed ./enh_exp_targetgene_exp | awk -v tissue="$i" 'BEGIN{OFS="\t"}{ if($13==tissue&&$12!="random"&&$14=="spe"){print $10} }' - |sort -u >./target_gene/TS_gene/$i"_TS_target"
  echo $i" 完成TS target gene"
done


#合并十个组织TS的靶基因到一个文件：
cd /home/ding/all_cmd/script/enh_statistics
awk 'BEGIN{OFS="\t";print "gene_name","tissue"}' >./target_gene/TS_gene/all_TS_target
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  awk -v tissue="$i" '{print $0,tissue}' ./target_gene/TS_gene/$i"_TS_target" >>./target_gene/TS_gene/all_TS_target
done













