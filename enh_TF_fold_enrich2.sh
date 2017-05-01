#!/bin/bash

#######################################################################
#TFBS_sorted.bed
#/home/ding/all_cmd/script/enh_statistics
cd /media/ding/000B49000006264C/eRNA_project/TFBS
#spe
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$7=="spe"){print $1,$2,$3,$4} }' /home/ding/all_cmd/script/enh_statistics/enh_spe  |closest-features --closest  --delim '\t' --dist - ./TFBS_sorted.bed |awk 'BEGIN{OFS="\t"}($NF<=2000&&$NF>=-2000){print $8}' - |sort |uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$1}' - >./a
#uni:
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$7=="uni"){print $1,$2,$3,$4} }' /home/ding/all_cmd/script/enh_statistics/enh_spe  |closest-features --closest  --delim '\t' --dist - ./TFBS_sorted.bed |awk 'BEGIN{OFS="\t"}($NF<=2000&&$NF>=-2000){print $8}' - |sort |uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$1}' - >./b
#other:
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$7=="other"){print $1,$2,$3,$4} }' /home/ding/all_cmd/script/enh_statistics/enh_spe  |closest-features --closest  --delim '\t' --dist - ./TFBS_sorted.bed |awk 'BEGIN{OFS="\t"}($NF<=2000&&$NF>=-2000){print $8}' - |sort |uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$1}' - >./c


#计算富集程度：
  awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$1]=$2;}  if(NR>FNR&&$1 in a){print $1,$2,a[$1]} if(NR>FNR&&!($1 in a)){print $1,$2,0.1} }' ./a ./c >./ac1
  awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$1]=$2;}  if(NR>FNR&&!($1 in a)){print $1,0.1,$2} }' ./c ./a >./ac2
  awk 'BEGIN{FS="\t";OFS="\t"}{print $0,$2/$3}' ./ac1 ./ac2 >./ac
  rm ./ac1 ./ac2
  
  awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$1]=$2;}  if(NR>FNR&&$1 in a){print $1,$2,a[$1]} if(NR>FNR&&!($1 in a)){print $1,$2,0.1} }' ./b ./c >./bc1
  awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$1]=$2;}  if(NR>FNR&&!($1 in a)){print $1,0.1,$2} }' ./c ./b >./bc2
  awk 'BEGIN{FS="\t";OFS="\t"}{print $0,$2/$3}' ./bc1 ./bc2  >./bc
  rm ./bc1 ./bc2
  
#然后bid/unbid/no_eRNA合并到一个文件并创建结果文件：
echo -e "TF\tspe_TF_fold_enrich\tuni_TF_fold_enrich" >/home/ding/all_cmd/script/enh_statistics/enh_TF_fold_enrich2
  awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1]=$NF} if(NR>FNR&&$1 in a){print $1,a[$1],$NF} }'  ./bc ./ac >>/home/ding/all_cmd/script/enh_statistics/enh_TF_fold_enrich2
  rm ./a ./b ./c ./ac ./bc
  echo "完成1"

#使用R做图：











#下面这个不减去other的背景：

#######################################################################
#TFBS_sorted.bed

#/home/ding/all_cmd/script/enh_statistics
cd /media/ding/000B49000006264C/eRNA_project/TFBS
#spe
closest-features --closest  --delim '\t' --dist  ./TFBS_sorted.bed /home/ding/all_cmd/script/enh_statistics/enh_spe_cluster.bed |awk 'BEGIN{OFS="\t"}($NF<=2000&&$NF>=-2000){print $4}' - |sort |uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$1}' - >./a
#uni:
closest-features --closest  --delim '\t' --dist ./TFBS_sorted.bed  /home/ding/all_cmd/script/enh_statistics/enh_uni_cluster.bed |awk 'BEGIN{OFS="\t"}($NF<=2000&&$NF>=-2000){print $4}' - |sort |uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$1}' - >./b
#other:
closest-features --closest  --delim '\t' --dist ./TFBS_sorted.bed /home/ding/all_cmd/script/enh_statistics/enh_other_cluster.bed  |awk 'BEGIN{OFS="\t"}($NF<=2000&&$NF>=-2000){print $4}' - |sort |uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$1}' - >./c


#直接将spe和uni的表达量/总数归一化后作为富集程度进行比较：
#spe:
awk 'BEGIN{FS="\t";OFS="\t";all=0}{all=all+$2;a[$1]=$2}END{ for(i in a){print i,a[i]/all*10^6} }' ./a  >./a_nor
#uni:
awk 'BEGIN{FS="\t";OFS="\t";all=0}{all=all+$2;a[$1]=$2}END{ for(i in a){print i,a[i]/all*10^6} }' ./b  >./b_nor


#然后eRNA/no_eRNA合并到一个文件并创建结果文件：
echo -e "TF\tspe_enrich\tuni_enrich" >/home/ding/all_cmd/script/enh_statistics/enh_TF_fold_enrich2
awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1]=$NF} if(NR>FNR&&$1 in a){print $1,a[$1],$NF} if(NR>FNR&&!($1 in a)){print $1,$2,0} }'  ./a_nor ./b_nor > ./ab_nor1
awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$1]=$2;}  if(NR>FNR&&!($1 in a)){print $1,0,$2} }'  ./b_nor ./a_nor > ./ab_nor2
awk 'BEGIN{FS="\t";OFS="\t"}{ print }'  ./ab_nor1  ./ab_nor2 >>/home/ding/all_cmd/script/enh_statistics/enh_TF_fold_enrich2
  
rm ./a ./b ./c  ./a_nor ./b_nor  ./ab_nor1 ./ab_nor2 
echo "完成1"














########################################################
########################################################
########################################################
#获得每个转录因子在每个组织中的增强子上的结合位点的数目作为其富集程度：

tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
for i in ${tissue[@]}
do 
  cd /home/ding/all_cmd/script/enh_statistics

   awk -v tissue=$i 'BEGIN{OFS="\t"}{ if(NR>1&&$8==tissue){printf "%s\t%d\t%d\t%s\n",$1,$2,$3,$4} }' ./enh_all_count_spe_melt|sortBed -i -| closest-features --closest  --delim '\t' --dist  /media/ding/000B49000006264C/eRNA_project/TFBS/TFBS_sorted.bed  - | awk 'BEGIN{OFS="\t"}($NF<=2000&&$NF>=-2000){print $4}' - |sort |uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$1}' - |sort -k 2 -nr |awk 'BEGIN{FS="\t";OFS="\t";all=0}{all=all+$2;a[$1]=$2}END{ for(i in a){printf "%s\t%g\n",i,a[i]} }' - >./TF_enrich/$i"_TF_enrich"

  echo $i"中的TF enrich完成！"
done


#不区分组织的转录因子富集程度：
cd /home/ding/all_cmd/script/enh_statistics
closest-features --closest  --delim '\t' --dist  /media/ding/000B49000006264C/eRNA_project/TFBS/TFBS_sorted.bed enh_all_cluster.bed | awk 'BEGIN{OFS="\t"}($NF<=2000&&$NF>=-2000){print $4}' - |sort |uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$1}' - |sort -k 2 -nr |awk 'BEGIN{FS="\t";OFS="\t";all=0}{all=all+$2;a[$1]=$2}END{ for(i in a){printf "%s\t%g\n",i,a[i]} }' - >./TF_enrich/all_TF_enrich


cd  ~/all_cmd/script/enh_statistics/TF_enrich
grep HNF4G ./* |awk 'BEGIN{OFS="\t"}{print $1,$2}' -|sort -n -k 2

#但是每个组织的转录因子富集程度得单独计算：



