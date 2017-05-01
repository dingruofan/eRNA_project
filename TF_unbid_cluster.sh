#!/bin/bash

cd /media/ding/000B49000006264C/eRNA_project/histone
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  #将各个组织唯一的TF保留到临时文件a
  awk '{print $10}' ./$i/enh_find_1/enh_unbid_TF.bed | sort -u >> /home/ding/all_cmd/script/enh_statistics/a
done
#gedit  /home/ding/all_cmd/script/enh_statistics/a

#将各个组织的TF去重复，保留唯一的TF到TF统计文件TF_unbid_count,并在首行添加列明：TF_all_tissue
awk '{print $0}' /home/ding/all_cmd/script/enh_statistics/a | sort -u |awk 'BEGIN{print "TF_all_tissue"}{print $0}' - > /home/ding/all_cmd/script/enh_statistics/TF_unbid_count
rm /home/ding/all_cmd/script/enh_statistics/a

#统计TF在各个组织中是否出现：
cd /media/ding/000B49000006264C/eRNA_project/histone
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  #将各个组织唯一的TF保留到临时文件a
  awk '{print $10}' ./$i/enh_find_1/enh_unbid_TF.bed | sort -u |awk 'BEGIN{} { if(NR==FNR){a[NR]=$1;next} if(NR>FNR){b[$1]} } END{ for(i=2;i<=length(a);i++){if(a[i] in b){print a[i]"\t"1}else{print  a[i]"\t"0}} }'  /home/ding/all_cmd/script/enh_statistics/TF_unbid_count - | awk -v tissue="$i" '{ if(NR==FNR){a[$1]=$2;next} if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }' - /home/ding/all_cmd/script/enh_statistics/TF_unbid_count > /home/ding/all_cmd/script/enh_statistics/a
  awk '{print $0}' /home/ding/all_cmd/script/enh_statistics/a > /home/ding/all_cmd/script/enh_statistics/TF_unbid_count
  rm /home/ding/all_cmd/script/enh_statistics/a
done

#统计每个TF在几个组织中存在
awk 'BEGIN{}{ if(NR>1){sum=0;for(i=2;i<=NF;i++){sum=$i+sum} print $1"\t"sum} }END{}' /home/ding/all_cmd/script/enh_statistics/TF_unbid_count |awk '{ if(NR==FNR){a[$1]=$2} if(NR>FNR&&FNR==1){print $0"\ttissue_statistics"} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }' - /home/ding/all_cmd/script/enh_statistics/TF_unbid_count >/home/ding/all_cmd/script/enh_statistics/a

awk '{print $0}' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/TF_unbid_count 
rm /home/ding/all_cmd/script/enh_statistics/a


cd /media/ding/000B49000006264C/eRNA_project/histone
#awk '{ if(NR==1){print} if($NF==10){print } }' /home/ding/all_cmd/script/enh_statistics/TF_unbid_count
awk '{ if(NR==1){print} if($NF==1){print } }' /home/ding/all_cmd/script/enh_statistics/TF_unbid_count

gedit /home/ding/all_cmd/script/enh_statistics/TF_unbid_count


