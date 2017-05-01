
#已经获得了/media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/cage_enh/bidir_pairs_expression_TPM.matrix双向每个转录本的表达量
#单向的：/media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/cage_enh/TCs_sorted_exp_TPM.matrix的转录本的表达量

#做bid eRNA表达量分析：
#从所有的enh_bid_cluster.bed匹配到每个组织的eRNA，然后匹配到表达矩阵。


#创建/home/ding/all_cmd/script/enh_statistics/enh_bid_cluster_for_exp.bed
  awk 'BEGIN{OFS="\t";}{ print $1,$2,$3,$4",-";print $1,$2,$3,$4",+";}' /home/ding/all_cmd/script/enh_statistics/enh_bid_cluster.bed >/home/ding/all_cmd/script/enh_statistics/enh_bid_cluster_for_exp.bed
  awk 'BEGIN{print "chr\tchr_start\tchr_end\tenhancer"}{print $0}' /home/ding/all_cmd/script/enh_statistics/enh_bid_cluster_for_exp.bed > /home/ding/all_cmd/script/enh_statistics/enh_bid_exp

tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  #首先获得每个组织的eRNA的表达量放在各个组织文件夹下，注意同一个enh的在不同组织的eRNA可能不同。
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  awk -v tissue=$i 'BEGIN{tissue_ID=0;OFS="\t"}{ if(NR==1){ for(i=1;i<=NF;i++){ if($i==tissue){tissue_ID=i} } } if((NR==FNR)&&(NR>1)){a[$1]=$tissue_ID}  if(NR>FNR&&($6",-" in a)){print $1,$2,$3,$4,$5,"-",$7,a[$6",-"]} if(NR>FNR&&($6",+" in a)){print $1,$2,$3,$4,$5,"+",$7,a[$6",+"]}}'  ../../fantom/All_tissue_CAGE/cage_enh/bidir_pairs_expression_TPM.matrix  ./enh_find/enh_bid.bed >./enh_find/enh_bid_exp.bed 
  
  #注意：上面这个双向的没有详细的eRNA的信息，下面是补充到另一个文件：./enh_find/enh_bid_exp2.bed 
cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  awk -v tissue=$i 'BEGIN{tissue_ID=0;OFS="\t"}{ if(NR==1){ for(i=1;i<=NF;i++){ if($i==tissue){tissue_ID=i} } } if((NR==FNR)&&(NR>1)){a[$1]=$tissue_ID}  if(NR>FNR&&($6",-" in a)){print $1,$2,$3,$4,"-",$6",-",$7,a[$6",-"]} if(NR>FNR&&($6",+" in a)){print $1,$2,$3,$4,"+",$6",+",$7,a[$6",+"]}}'  ../../fantom/All_tissue_CAGE/cage_enh/bidir_pairs_expression_TPM.matrix  ./enh_find/enh_bid.bed >./enh_find/enh_bid_exp2.bed 

 #然后再将表达量合并到enh_bid_exp，没有的表达量为0。
  awk -v tissue=$i '{   if(NR==FNR){ a[$4","$6]=$NF }  if(NR>FNR){ if($4 in a){ print  $1,$2,$3,$4,a[$4] }else{print  $1,$2,$3,$4,0} }   }' ./enh_find/enh_bid_exp.bed   /home/ding/all_cmd/script/enh_statistics/enh_bid_cluster_for_exp.bed |awk  -v tissue=$i 'BEGIN{OFS="\t"}{  if(NR==FNR){a[$1"\t"$2"\t"$3"\t"$4]=$5} if(NR>FNR&&FNR==1){print $0"\t"tissue}  if(NR>FNR&&FNR>1){print $0,a[$1"\t"$2"\t"$3"\t"$4]} }' -  /home/ding/all_cmd/script/enh_statistics/enh_bid_exp >/home/ding/all_cmd/script/enh_statistics/a 
  
  awk '{ print $0 }' /home/ding/all_cmd/script/enh_statistics/a > /home/ding/all_cmd/script/enh_statistics/enh_bid_exp
  rm /home/ding/all_cmd/script/enh_statistics/a
  echo $i"完成1"
done

#————————————————————————————————————————————————————————————————————————————————————————————————————————————————————


#做unbid eRAN的表达量分析：
#从所有的enh_unbid_cluster.bed匹配到每个组织的eRNA，然后匹配到表达矩阵。

awk 'BEGIN{print "chr\tchr_start\tchr_end\tenhancer"}{print $0}' /home/ding/all_cmd/script/enh_statistics/enh_unbid_cluster.bed > /home/ding/all_cmd/script/enh_statistics/enh_unbid_exp

tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  #首先获得每个组织的eRNA的表大两放在各个组织文件夹下，注意同一个enh的在不同组织的eRNA可能不同。
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  awk -v tissue=$i 'BEGIN{tissue_ID=0 }{ if(NR==1){ for(i=1;i<=NF;i++){ if($i==tissue){tissue_ID=i} } } if((NR==FNR)&&(NR>1)){a[$1]=$tissue_ID}  if(NR>FNR){print $0"\t"a[$6]} }'  ../../fantom/All_tissue_CAGE/cage_enh/TCs_sorted_exp_TPM.matrix  ./enh_find_1/enh_unbid.bed  >./enh_find_1/enh_unbid_exp.bed 

 #然后再将表达量合并到enh_bid_exp，没有的表达量为0。
  awk -v tissue=$i '{   if(NR==FNR){ a[$1"\t"$2"\t"$3"\t"$4]=$8 }  if(NR>FNR){ if($1"\t"$2"\t"$3"\t"$4 in a){ print  $1"\t"$2"\t"$3"\t"$4"\t"a[$1"\t"$2"\t"$3"\t"$4] }else{print  $1"\t"$2"\t"$3"\t"$4"\t0"} }   }' ./enh_find_1/enh_unbid_exp.bed   /home/ding/all_cmd/script/enh_statistics/enh_unbid_cluster.bed |awk  -v tissue=$i '{  if(NR==FNR){a[$1"\t"$2"\t"$3"\t"$4]=$5} if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1"\t"$2"\t"$3"\t"$4]}  }' -  /home/ding/all_cmd/script/enh_statistics/enh_unbid_exp >/home/ding/all_cmd/script/enh_statistics/a
  awk '{ print $0 }' /home/ding/all_cmd/script/enh_statistics/a > /home/ding/all_cmd/script/enh_statistics/enh_unbid_exp 
  rm /home/ding/all_cmd/script/enh_statistics/a

  echo $i"完成2"
done














