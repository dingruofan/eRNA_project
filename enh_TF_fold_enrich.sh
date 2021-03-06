#!/bin/bash

#######################################################################
# 
# #先计算bid 每个TF的TFBS数目,富集程度为TFBS/sum(TFBS)
# tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
# #tissue=(Adrenal)
# for i in ${tissue[@]}
# do 
# cd /media/ding/000B49000006264C/eRNA_project/histone/$i
# #对bid/unbid/no_eRNA的TF做富集：
# awk 'BEGIN{FS="\t";OFS="\t"}{ print $10 }' ./enh_find/enh_bid_TF.bed | sort |uniq -c |awk 'BEGIN{OFS="\t";sum=0;}{ sum+=$1;a[NR]=$1;b[NR]=$2"\t"$1 }END{ for(i=1;i<=length(a);i++){print b[i],a[i]/sum*1000} }' -  >./a
# awk 'BEGIN{FS="\t";OFS="\t"}{ print $10 }' ./enh_find_1/enh_unbid_TF.bed | sort |uniq -c |awk 'BEGIN{OFS="\t";sum=0;}{ sum+=$1;a[NR]=$1;b[NR]=$2"\t"$1 }END{ for(i=1;i<=length(a);i++){print b[i],a[i]/sum*1000} }' -  >./b
# awk 'BEGIN{FS="\t";OFS="\t"}{ print $8 }' ./enh_find/enh_no_eRNA_named_TFBS_exp | sort |uniq -c |awk 'BEGIN{OFS="\t";sum=0;}{ sum+=$1;a[NR]=$1;b[NR]=$2"\t"$1 }END{ for(i=1;i<=length(a);i++){print b[i],a[i]/sum*1000} }' -  >./c
# #合并bid和no_eRNA；以及unbid和no_eRNA，并计算fold enrichment:
# awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$1]=$2"\t"$3;}  if(NR>FNR&&$1 in a){print $0,a[$1]} }' ./a ./c | awk '{print $0,$3-$5}' - >./ac
# awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$1]=$2"\t"$3;}  if(NR>FNR&&$1 in a){print $0,a[$1]} }' ./b ./c | awk '{print $0,$3-$5}' - >./bc
# 
# #然后每个组织的bid/unbid/no_eRNA合并到一个文件：
# awk -v tissue=$i '{ if(NR==FNR){a[$1]=$NF} if(NR>FNR&&$1 in a){print $1,a[$1],$NF,tissue} }' ./ac ./bc >./enh_TF_fold_enrich
# rm ./a ./b ./c ./ac ./bc
# 
# echo $i"完成1"
# done
# 
# 
# #创建结果文件：
# echo -e "TF\tbid_TF_fold_enrich\tunbid_TF_fold_enrich\ttissue" >/home/ding/all_cmd/script/enh_statistics/enh_TF_fold_enrich
# 
# #合并各个组织的/enh_TF_fold_enrich到一个文件：
# tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
# #tissue=(Adrenal)
# for i in ${tissue[@]}
# do 
# cd /media/ding/000B49000006264C/eRNA_project/histone/$i
# awk '{print $0}' ./enh_TF_fold_enrich >>/home/ding/all_cmd/script/enh_statistics/enh_TF_fold_enrich
# echo $i"完成2"
# done
# 
# #使用R做图：


################################################
#上面的有点问题：不应该按照每个组织的TF进行统计，因为TF在所有组织中是一样的
#TFBS_sorted.bed
#/home/ding/all_cmd/script/enh_statistics
cd /media/ding/000B49000006264C/eRNA_project/TFBS
#bid
closest-features --closest  --delim '\t' --dist /home/ding/all_cmd/script/enh_statistics/enh_bid_cluster.bed ./TFBS_sorted.bed |awk 'BEGIN{OFS="\t"}($NF<=2000&&$NF>=-2000){print $8}' - |sort |uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$1}' - >./a
#unbid:
closest-features --closest  --delim '\t' --dist /home/ding/all_cmd/script/enh_statistics/enh_unbid_cluster.bed ./TFBS_sorted.bed |awk 'BEGIN{OFS="\t"}($NF<=2000&&$NF>=-2000){print $8}' - |sort |uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$1}' - >./b
#no_eRNA:
closest-features --closest  --delim '\t' --dist /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_cluster.bed ./TFBS_sorted.bed |awk 'BEGIN{OFS="\t"}($NF<=2000&&$NF>=-2000){print $8}' - |sort |uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$1}' - >./c


#计算富集程度：
  awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$1]=$2;}  if(NR>FNR&&$1 in a){print $0,a[$1]} }' ./a ./c | awk '{print $0"\t"$2/$3}' - >./ac
  awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$1]=$2;}  if(NR>FNR&&$1 in a){print $0,a[$1]} }' ./b ./c | awk '{print $0"\t"$2/$3}' - >./bc
  
#然后bid/unbid/no_eRNA合并到一个文件并创建结果文件：
echo -e "TF\tbid_TF_fold_enrich\tunbid_TF_fold_enrich" >/home/ding/all_cmd/script/enh_statistics/enh_TF_fold_enrich
  awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1]=$NF} if(NR>FNR&&$1 in a){print $1,a[$1],$NF} }'  ./bc ./ac >>/home/ding/all_cmd/script/enh_statistics/enh_TF_fold_enrich
  rm ./a ./b ./c ./ac ./bc
  echo "完成1"

#使用R做图：



###################################################################
###################################################################
#不与背景比较：
#/home/ding/all_cmd/script/enh_statistics
cd /media/ding/000B49000006264C/eRNA_project/TFBS
#bid
closest-features --closest  --delim '\t' --dist  ./TFBS_sorted.bed  /home/ding/all_cmd/script/enh_statistics/enh_bid_cluster.bed |awk 'BEGIN{OFS="\t"}($NF<=2000&&$NF>=-2000){print $4}' - |sort |uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$1}' - >./a
#unbid:
closest-features --closest  --delim '\t' --dist  ./TFBS_sorted.bed  /home/ding/all_cmd/script/enh_statistics/enh_unbid_cluster.bed|awk 'BEGIN{OFS="\t"}($NF<=2000&&$NF>=-2000){print $4}' - |sort |uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$1}' - >./b
#no_eRNA:
closest-features --closest  --delim '\t' --dist  ./TFBS_sorted.bed  /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_cluster.bed  |awk 'BEGIN{OFS="\t"}($NF<=2000&&$NF>=-2000){print $4}' - |sort |uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$1}' - >./c

#合并bid和unbid的作为eRNAs-enh，并使用TFBS/总数归一化
awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$1]=$2;}  if(NR>FNR&&$1 in a){print $1,$2,a[$1]} if(NR>FNR&&!($1 in a)){print $1,$2,0} }' ./a ./b >./ab1
awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$1]=$2;}  if(NR>FNR&&!($1 in a)){print $1,0,$2} }' ./b ./a >./ab2
awk 'BEGIN{FS="\t";OFS="\t";all=0;}{ all=all+$2+$3; a[$1]=$2+$3; }END{ for(i in a){ print i,a[i]/all*10^6 } }' ./ab1 ./ab2 >./ab_nor

#no-eRNAs归一化：
awk 'BEGIN{FS="\t";OFS="\t";all=0;}{ all=all+$2;a[$1]=$2; }END{ for(i in a){ print i,a[i]/all*10^6 } } '  ./c > ./c_nor



#然后eRNA/no_eRNA合并到一个文件并创建结果文件：
echo -e "TF\teRNA_enrich\tno_eRNA_enrich" >/home/ding/all_cmd/script/enh_statistics/enh_TF_fold_enrich
  awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1]=$NF} if(NR>FNR&&$1 in a){print $1,a[$1],$NF} if(NR>FNR&&!($1 in a)){print $1,$2,0} }'  ./ab_nor ./c_nor > ./abc_nor1
  awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$1]=$2;}  if(NR>FNR&&!($1 in a)){print $1,0,$2} }'  ./c_nor ./ab_nor > ./abc_nor2
  awk 'BEGIN{FS="\t";OFS="\t"}{ print }'  ./abc_nor1  ./abc_nor2 >>/home/ding/all_cmd/script/enh_statistics/enh_TF_fold_enrich
  
  rm ./a ./b ./c  ./ab1 ./ab2 ./ab_nor ./c_nor ./abc_nor1 ./abc_nor2 
  echo "完成1"













