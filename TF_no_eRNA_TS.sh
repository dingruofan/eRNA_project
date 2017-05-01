#!/bin/bash
cd  /media/ding/000B49000006264C/eRNA_project/TF_exp
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)

#建立要输出的文件，并输入固定列TF：
awk '{print $1}' /home/ding/all_cmd/script/enh_statistics/TF_no_eRNA_count > /home/ding/all_cmd/script/enh_statistics/TF_no_eRNA_count_exp

#根据每个组织的转录因子表大量文件，循环输入转录因子的表大量：
#如果想将表达量为0的替换为相应组织的表达量，将{print $1"\t"$tissue_id}替换为 ：{print $1"\t"a[$1]}
for i in ${tissue[@]}
do 
  
  awk -v tissue="$i" 'BEGIN{tissue_id=0} {  if(NR==FNR){a[$2]=$3;next} if(NR>FNR&&FNR==1){for(j=1;j<=NF;j++){if($j==tissue){tissue_id=j}} }   if(NR>FNR&&FNR>1&&$tissue_id==0){print $1"\t"a[$1]} if(NR>FNR&&FNR>1&&$tissue_id!=0){print $1"\t"a[$1]} } ' ./$i/TF_exp  /home/ding/all_cmd/script/enh_statistics/TF_no_eRNA_count |awk -v tissue="$i" '{ if(NR==FNR){a[$1]=$2;next} if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }' - /home/ding/all_cmd/script/enh_statistics/TF_no_eRNA_count_exp > /home/ding/all_cmd/script/enh_statistics/a
  awk '{print $0}' /home/ding/all_cmd/script/enh_statistics/a >  /home/ding/all_cmd/script/enh_statistics/TF_no_eRNA_count_exp
  rm  /home/ding/all_cmd/script/enh_statistics/a

done



#去掉有NA的行：

awk '{ if(NR==1){print $0} if(NR>1&&!match($0,"\tNA")){print $0} }' /home/ding/all_cmd/script/enh_statistics/TF_no_eRNA_count_exp > /home/ding/all_cmd/script/enh_statistics/a
awk '{print $0}' /home/ding/all_cmd/script/enh_statistics/a >  /home/ding/all_cmd/script/enh_statistics/TF_no_eRNA_count_exp
rm  /home/ding/all_cmd/script/enh_statistics/a

#调用R脚本计算每个TF的特异性分数并添加到最后一列，还是写入原文件：
Rscript /home/ding/all_cmd/script/TF_no_eRNA_tissue_score.R


gedit /home/ding/all_cmd/script/enh_statistics/TF_no_eRNA_count_exp


