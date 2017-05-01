#!/bin/bash

#____________________________________________________________________________
#eRNA表达量与转录因子结合为点个数分析：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
for i in ${tissue[@]}
do 
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  #获得每个组织的enh和他所对应的转录因子TFBS个数。
  #bid:
  awk -v tissue=$i '{  if(NR==FNR){a[$4]++;next}  if(NR>FNR){ if($4 in a){print $4"\t"$8"\t"a[$4]"\t"tissue}else{print $4"\t"$8"\t0\t"tissue} } }' ./enh_find/enh_bid_TF.bed   ./enh_find/enh_bid_exp.bed >./enh_find/enh_bid_exp_TF
  #unbid:
  awk -v tissue=$i '{  if(NR==FNR){a[$4]++;next}  if(NR>FNR){ if($4 in a){print $4"\t"$8"\t"a[$4]"\t"tissue}else{print $4"\t"$8"\t0\t"tissue} } }' ./enh_find_1/enh_unbid_TF.bed   ./enh_find_1/enh_unbid_exp.bed >./enh_find_1/enh_unbid_exp_TF
  #no_eRNA:没有表达量。。。
  awk -v tissue=$i '{  if(NR==FNR){a[$4]++;next}  if(NR>FNR){ if($4 in a){print $4"\t0""\t"a[$4]"\t"tissue}else{print $4"\t0""\t0\t"tissue} } }' ./enh_find/enh_no_eRNA_named_TFBS_exp    ./enh_find/enh_no_eRNA_named.bed >./enh_find/enh_no_eRNA_exp_TF

  echo $i"完成1"
done

###后面补充的，对下面的无影响：
########################################
#获得10个组织的所有转录因子结合位点及enh的所有数据：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  #获得每个组织的enh和他所对应的所有转录因子统计：
  #包括bid unbid 和no-eRNA的
  awk -v tissue=$i 'BEGIN{OFS="\t"}{  if(FILENAME=="./enh_find/enh_bid_TF.bed"){print $1,$2,$3,$4,$7,$8,$9,$10}  if(FILENAME=="./enh_find_1/enh_unbid_TF.bed"){print $1,$2,$3,$4,$7,$8,$9,$10} if(FILENAME=="./enh_find/enh_no_eRNA_named_TFBS_exp"){print $1,$2,$3,$4,$5,$6,$7,$8}  }' ./enh_find/enh_bid_TF.bed ./enh_find_1/enh_unbid_TF.bed  ./enh_find/enh_no_eRNA_named_TFBS_exp >/home/ding/all_cmd/script/enh_statistics/enh_TF/"enh_TF_"$i

  echo $i"完成1"
done

#然后将上面的10个组织的enh-TF数据合并为一个文件：
awk 'BEGIN{OFS="\t";print "enh_chr","enh_start","enh_end","enh_name","TF_chr","TF_start","TF_end","TF_name","tissue"}' >/home/ding/all_cmd/script/enh_statistics/enh_TF/enh_TF_all
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  awk -v tissue=$i 'BEGIN{OFS="\t"} {print $0,tissue}' /home/ding/all_cmd/script/enh_statistics/enh_TF/"enh_TF_"$i >>/home/ding/all_cmd/script/enh_statistics/enh_TF/enh_TF_all
done

#使用R进行整理：做成类似spe_TF_enh_tissue_3这样的格式：
#使用enh_ts3.R整理：
Rscript  enh_ts3.R
########################################


#____________________________________________________________________________
#bid/unbid/no_eRNA的exp_TF分别合并到all_cmd里：
cd /home/ding/all_cmd/script 
 awk 'BEGIN{print "enhancer\teRNA_exp_TPM\tTFBS_num\ttissue"}' > /home/ding/all_cmd/script/enh_statistics/enh_bid_exp_TF_all  

tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
for i in ${tissue[@]}
do 
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  awk '{  if(NR==FNR){print $0} if(NR>FNR){print $0} }'  /home/ding/all_cmd/script/enh_statistics/enh_bid_exp_TF_all  ./enh_find/enh_bid_exp_TF >  /home/ding/all_cmd/script/enh_statistics/a
  awk '{  print $0  }' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_bid_exp_TF_all 
  rm /home/ding/all_cmd/script/enh_statistics/a

  echo $i"完成2"
done
#____________________________________________________________________________
#unbid:
cd /home/ding/all_cmd/script 
 awk 'BEGIN{print "enhancer\teRNA_exp_TPM\tTFBS_num\ttissue"}' > /home/ding/all_cmd/script/enh_statistics/enh_unbid_exp_TF_all  

tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
for i in ${tissue[@]}
do 
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  awk '{  if(NR==FNR){print $0} if(NR>FNR){print $0} }'  /home/ding/all_cmd/script/enh_statistics/enh_unbid_exp_TF_all  ./enh_find_1/enh_unbid_exp_TF >  /home/ding/all_cmd/script/enh_statistics/a
  awk '{  print $0  }' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_unbid_exp_TF_all 
  rm /home/ding/all_cmd/script/enh_statistics/a

  echo $i"完成3"
done
#____________________________________________________________________________
#no_eRNA:
cd /home/ding/all_cmd/script 
 awk 'BEGIN{print "enhancer\teRNA_exp_TPM\tTFBS_num\ttissue"}' > /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_exp_TF_all  

tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
for i in ${tissue[@]}
do 
 
  awk '{  if(NR==FNR){print $0} if(NR>FNR){print $0} }'  /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_exp_TF_all  ./enh_find/enh_no_eRNA_exp_TF >  /home/ding/all_cmd/script/enh_statistics/a
  awk '{  print $0  }' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_exp_TF_all 
  rm /home/ding/all_cmd/script/enh_statistics/a

  echo $i"完成4"
done

#____________________________________________________________________________
#每个组织的random位点的TFBS结合位点：
#random位点的个数与这个组织增强子的数目一致：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
for i in ${tissue[@]}
do
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  #获得这个组织增强子的个数：result.bed
  enhancer_num=`wc -l ./enh_find/result.bed|awk '{ print $1 }' -`
  bedtools random -l 500 -n $enhancer_num -seed 2 -g  /home/ding/all_cmd/hg19.chrom_24.sizes |awk '{  print $1"\t"$2"\t"$3"\trandom_"$4  }' - |bedops --not-element-of - /media/ding/000B49000006264C/eRNA_project/gencode/protein_gene.bed |sortBed -i - >./enh_find/randomBed_cluster.bed
  #获得2kb附近的转录因子：
  closest-features --closest  --delim '\t' --dist ../../TFBS/TFBS_sorted.bed ./enh_find/randomBed_cluster.bed |awk '($NF<=2000&&$NF>=-2000){print }' - |awk -v tissue=$i 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$9]++} if(NR>FNR&&$4 in a){print $4,0,a[$4],tissue} if((NR>FNR)&&!($4 in a)){print $4,0,0,tissue} }' -  ./enh_find/randomBed_cluster.bed >./enh_find/enh_random_exp_TF
  
  echo $i"完成5"
done

#random
#将各个组织random的分别合并到all_cmd里：
cd /home/ding/all_cmd/script 
 awk 'BEGIN{print "enhancer\teRNA_exp_TPM\tTFBS_num\ttissue"}' > /home/ding/all_cmd/script/enh_statistics/enh_random_exp_TF_all  

tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
for i in ${tissue[@]}
do 
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  awk '{  if(NR==FNR){print $0} if(NR>FNR){print $0} }'  /home/ding/all_cmd/script/enh_statistics/enh_random_exp_TF_all  ./enh_find/enh_random_exp_TF >  /home/ding/all_cmd/script/enh_statistics/a
  awk '{  print $0  }' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_random_exp_TF_all 
  rm /home/ding/all_cmd/script/enh_statistics/a

  echo $i"完成6"
done










