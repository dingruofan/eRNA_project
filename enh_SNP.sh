#!/bin/bash

#######################################################################
#处理下载的SNP147：


#产生随机位点进行对照比较：
#我的enh的平均长度约为400bp。
cd /home/ding/all_cmd/script
bedtools random -l 400 -n 50000 -seed 2 -g  ../hg19.chrom_24.sizes |awk '{  print $1"\t"$2"\t"$3"\trandom_"$4  }' - |bedops --not-element-of - /media/ding/000B49000006264C/eRNA_project/gencode/protein_coding_gene_up1000 |sortBed -i - >./enh_statistics/randomBed_cluster.bed


#比较bid/unbid/no_eRNA 2kb附近的SNP结合位点的多少。
#并注释特异性的enh:通过 /home/ding/all_cmd/script/enh_statistics/enh_spe
#####特异性必须先运行
 /home/ding/all_cmd/script/enh_ts.R

cd /media/ding/000B49000006264C/eRNA_project

#保留以下列："name","strand","observed","class","func","enh","dist","type","isspe"
#注意：去掉了gene外的SNP位点：awk 'NR>1{print $1"\t"$3"\t"$4}' /media/ding/000B49000006264C/eRNA_project/gencode/protein_gene 
#bid:
  closest-features --closest  --delim '\t' --dist ./SNP/SNP147_sorted.bed  /home/ding/all_cmd/script/enh_statistics/enh_bid_cluster.bed | bedops --not-element-of 1 - /media/ding/000B49000006264C/eRNA_project/gencode/protein_gene.bed |awk 'BEGIN{OFS="\t"}($NF<=0&&$NF>=0){print $4,$6,$7,$8,$9,$13,$NF,"bid"}' - |awk 'BEGIN{OFS="\t"}{ if(NR==FNR&&NR>1){a[$4]=$NF}  if(NR>FNR){print $0,a[$6]} }' /home/ding/all_cmd/script/enh_statistics/enh_spe - > /home/ding/all_cmd/script/enh_statistics/enh_bid_cluster_SNP.bed
#unbid:
  closest-features --closest  --delim '\t' --dist ./SNP/SNP147_sorted.bed  /home/ding/all_cmd/script/enh_statistics/enh_unbid_cluster.bed | bedops --not-element-of 1 - /media/ding/000B49000006264C/eRNA_project/gencode/protein_gene.bed |awk 'BEGIN{OFS="\t"}($NF<=0&&$NF>=0){print $4,$6,$7,$8,$9,$13,$NF,"unbid"}' - |awk 'BEGIN{OFS="\t"}{ if(NR==FNR&&NR>1){a[$4]=$NF}  if(NR>FNR){print $0,a[$6]} }' /home/ding/all_cmd/script/enh_statistics/enh_spe - > /home/ding/all_cmd/script/enh_statistics/enh_unbid_cluster_SNP.bed
#no_eRNA:
 closest-features --closest  --delim '\t' --dist ./SNP/SNP147_sorted.bed  /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_cluster.bed | bedops --not-element-of 1 - /media/ding/000B49000006264C/eRNA_project/gencode/protein_gene.bed  |awk 'BEGIN{OFS="\t"}($NF<=0&&$NF>=0){print $4,$6,$7,$8,$9,$13,$NF,"no_eRNA"}' -  |awk 'BEGIN{OFS="\t"}{ if(NR==FNR&&NR>1){a[$4]=$NF}  if(NR>FNR){print $0,a[$6]} }' /home/ding/all_cmd/script/enh_statistics/enh_spe - > /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_cluster_SNP.bed
#randomBed:
 closest-features --closest  --delim '\t' --dist ./SNP/SNP147_sorted.bed  /home/ding/all_cmd/script/enh_statistics/randomBed_cluster.bed | bedops --not-element-of 1 - /media/ding/000B49000006264C/eRNA_project/gencode/protein_gene.bed  |awk 'BEGIN{OFS="\t"}($NF<=0&&$NF>=0){print $4,$6,$7,$8,$9,$13,$NF,"random"}' - > /home/ding/all_cmd/script/enh_statistics/randomBed_cluster_SNP.bed

#统计bid/unbid/no_eRNA/random的每个增强子的SNP位点数目：
#bid:
cd /home/ding/all_cmd/script/enh_statistics
  awk '{  if(NR==FNR){ a[$6]++;next }   if(NR>FNR){ if($4 in a){ print $4"\t"a[$4] }else{print $4"\t"0} } }'  ./enh_bid_cluster_SNP.bed  ./enh_bid_cluster.bed >./enh_bid_cluseter_SNP_num
#unbid:
  awk '{  if(NR==FNR){ a[$6]++;next }   if(NR>FNR){ if($4 in a){ print $4"\t"a[$4] }else{print $4"\t"0} } }'  ./enh_unbid_cluster_SNP.bed  ./enh_unbid_cluster.bed >./enh_unbid_cluseter_SNP_num
#no_eRNA
  awk '{  if(NR==FNR){ a[$6]++;next }   if(NR>FNR){ if($4 in a){ print $4"\t"a[$4] }else{print $4"\t"0} } }'  ./enh_no_eRNA_cluster_SNP.bed  ./enh_no_eRNA_cluster.bed >./enh_no_eRNA_cluseter_SNP_num
#randomBed
  awk '{  if(NR==FNR){ a[$6]++;next }   if(NR>FNR){ if($4 in a){ print $4"\t"a[$4] }else{print $4"\t"0} } }'  ./randomBed_cluster_SNP.bed  ./randomBed_cluster.bed >./randomBed_cluster_SNP_num



#后期补充：
#文件格式是：enh single insertion deletion
cd /home/ding/all_cmd/script/enh_statistics
#bid
awk 'BEGIN{OFS="\t"}{  if(NR==FNR&&$4=="single"){a[$6]++;next} if(NR==FNR&&$4=="insertion"){b[$6]++;next} if(NR==FNR&&$4=="deletion"){c[$6]++;next}  if(NR>FNR){ if(a[$4]==""){a[$4]=0} if(b[$4]==""){b[$4]=0} if(c[$4]==""){c[$4]=0} ;print $4,a[$4],b[$4],c[$4] } }'  ./enh_bid_cluster_SNP.bed  ./enh_bid_cluster.bed >./enh_bid_cluseter_SNP_num_detail
#unbid
awk 'BEGIN{OFS="\t"}{  if(NR==FNR&&$4=="single"){a[$6]++;next} if(NR==FNR&&$4=="insertion"){b[$6]++;next} if(NR==FNR&&$4=="deletion"){c[$6]++;next}  if(NR>FNR){ if(a[$4]==""){a[$4]=0} if(b[$4]==""){b[$4]=0} if(c[$4]==""){c[$4]=0} ;print $4,a[$4],b[$4],c[$4] } }'  ./enh_unbid_cluster_SNP.bed  ./enh_unbid_cluster.bed >./enh_unbid_cluseter_SNP_num_detail
#no_eRNA
awk 'BEGIN{OFS="\t"}{  if(NR==FNR&&$4=="single"){a[$6]++;next} if(NR==FNR&&$4=="insertion"){b[$6]++;next} if(NR==FNR&&$4=="deletion"){c[$6]++;next}  if(NR>FNR){ if(a[$4]==""){a[$4]=0} if(b[$4]==""){b[$4]=0} if(c[$4]==""){c[$4]=0} ;print $4,a[$4],b[$4],c[$4] } }'  ./enh_no_eRNA_cluster_SNP.bed  ./enh_no_eRNA_cluster.bed >./enh_no_eRNA_cluseter_SNP_num_detail
#randomBed
awk 'BEGIN{OFS="\t"}{  if(NR==FNR&&$4=="single"){a[$6]++;next} if(NR==FNR&&$4=="insertion"){b[$6]++;next} if(NR==FNR&&$4=="deletion"){c[$6]++;next}  if(NR>FNR){ if(a[$4]==""){a[$4]=0} if(b[$4]==""){b[$4]=0} if(c[$4]==""){c[$4]=0} ;print $4,a[$4],b[$4],c[$4] } }'  ./randomBed_cluster_SNP.bed  ./randomBed_cluster.bed >./randomBed_cluster_SNP_num_detail

#后期补充：
#spe:
cd /home/ding/all_cmd/script/enh_statistics
  awk 'BEGIN{OFS="\t"}{  if(NR==FNR&&$4=="single"){a[$6]++;next} if(NR==FNR&&$4=="insertion"){b[$6]++;next} if(NR==FNR&&$4=="deletion"){c[$6]++;next}  if(NR>FNR){ if(a[$4]==""){a[$4]=0} if(b[$4]==""){b[$4]=0} if(c[$4]==""){c[$4]=0} ;print $4,a[$4],b[$4],c[$4] } }' ./enh_spe_cluster_SNP.bed  ./enh_spe_cluster.bed >./enh_spe_cluseter_SNP_num_detail
#uni:
  awk 'BEGIN{OFS="\t"}{  if(NR==FNR&&$4=="single"){a[$6]++;next} if(NR==FNR&&$4=="insertion"){b[$6]++;next} if(NR==FNR&&$4=="deletion"){c[$6]++;next}  if(NR>FNR){ if(a[$4]==""){a[$4]=0} if(b[$4]==""){b[$4]=0} if(c[$4]==""){c[$4]=0} ;print $4,a[$4],b[$4],c[$4] } }' ./enh_uni_cluster_SNP.bed  ./enh_uni_cluster.bed >./enh_uni_cluseter_SNP_num_detail
#other
 awk 'BEGIN{OFS="\t"}{  if(NR==FNR&&$4=="single"){a[$6]++;next} if(NR==FNR&&$4=="insertion"){b[$6]++;next} if(NR==FNR&&$4=="deletion"){c[$6]++;next}  if(NR>FNR){ if(a[$4]==""){a[$4]=0} if(b[$4]==""){b[$4]=0} if(c[$4]==""){c[$4]=0} ;print $4,a[$4],b[$4],c[$4] } }' ./enh_other_cluster_SNP.bed  ./enh_other_cluster.bed >./enh_other_cluseter_SNP_num_detail
#random
  awk 'BEGIN{OFS="\t"}{  if(NR==FNR&&$4=="single"){a[$6]++;next} if(NR==FNR&&$4=="insertion"){b[$6]++;next} if(NR==FNR&&$4=="deletion"){c[$6]++;next}  if(NR>FNR){ if(a[$4]==""){a[$4]=0} if(b[$4]==""){b[$4]=0} if(c[$4]==""){c[$4]=0} ;print $4,a[$4],b[$4],c[$4] } }' ./randomBed_cluster_SNP.bed  ./randomBed_cluster.bed >./randomBed_cluster_SNP_num_detail









###################################################################
#统计bid/unbid/no_eRNA/random的
cd /media/ding/000B49000006264C/eRNA_project
#妈的，isspe类型的bed文件忘记排序了。。
#spe:
sortBed -i  /home/ding/all_cmd/script/enh_statistics/enh_spe_cluster.bed | closest-features --closest  --delim '\t' --dist ./SNP/SNP147_sorted.bed - | bedops --not-element-of 1 - /media/ding/000B49000006264C/eRNA_project/gencode/protein_gene.bed |awk 'BEGIN{OFS="\t"}($NF<=0&&$NF>=0){print $4,$6,$7,$8,$9,$13,$NF,"spe"}' - |awk 'BEGIN{OFS="\t"}{ if(NR==FNR&&NR>1){a[$4]=$5}  if(NR>FNR){print $0,a[$6]} }' /home/ding/all_cmd/script/enh_statistics/enh_spe - > /home/ding/all_cmd/script/enh_statistics/enh_spe_cluster_SNP.bed
#uni
sortBed -i /home/ding/all_cmd/script/enh_statistics/enh_uni_cluster.bed | closest-features --closest  --delim '\t' --dist ./SNP/SNP147_sorted.bed  - | bedops --not-element-of 1 - /media/ding/000B49000006264C/eRNA_project/gencode/protein_gene.bed |awk 'BEGIN{OFS="\t"}($NF<=0&&$NF>=0){print $4,$6,$7,$8,$9,$13,$NF,"uni"}' - |awk 'BEGIN{OFS="\t"}{ if(NR==FNR&&NR>1){a[$4]=$5}  if(NR>FNR){print $0,a[$6]} }' /home/ding/all_cmd/script/enh_statistics/enh_spe - > /home/ding/all_cmd/script/enh_statistics/enh_uni_cluster_SNP.bed
#other:
sortBed -i /home/ding/all_cmd/script/enh_statistics/enh_other_cluster.bed | closest-features --closest  --delim '\t' --dist ./SNP/SNP147_sorted.bed  - | bedops --not-element-of 1 - /media/ding/000B49000006264C/eRNA_project/gencode/protein_gene.bed |awk 'BEGIN{OFS="\t"}($NF<=0&&$NF>=0){print $4,$6,$7,$8,$9,$13,$NF,"other"}' - |awk 'BEGIN{OFS="\t"}{ if(NR==FNR&&NR>1){a[$4]=$5}  if(NR>FNR){print $0,a[$6]} }' /home/ding/all_cmd/script/enh_statistics/enh_spe - > /home/ding/all_cmd/script/enh_statistics/enh_other_cluster_SNP.bed
#randomBed:
 closest-features --closest  --delim '\t' --dist ./SNP/SNP147_sorted.bed  /home/ding/all_cmd/script/enh_statistics/randomBed_cluster.bed | bedops --not-element-of 1 - /media/ding/000B49000006264C/eRNA_project/gencode/protein_gene.bed  |awk 'BEGIN{OFS="\t"}($NF<=0&&$NF>=0){print $4,$6,$7,$8,$9,$13,$NF,"random"}' - > /home/ding/all_cmd/script/enh_statistics/randomBed_cluster_SNP.bed




#统计spe/uni/other/random的每个增强子的SNP位点数目：
#spe:
cd /home/ding/all_cmd/script/enh_statistics
  awk '{  if(NR==FNR){ a[$6]++;next }   if(NR>FNR){ if($4 in a){ print $4"\t"a[$4] }else{print $4"\t"0} } }'  ./enh_spe_cluster_SNP.bed  ./enh_spe_cluster.bed >./enh_spe_cluseter_SNP_num
#uni:
  awk '{  if(NR==FNR){ a[$6]++;next }   if(NR>FNR){ if($4 in a){ print $4"\t"a[$4] }else{print $4"\t"0} } }'  ./enh_uni_cluster_SNP.bed  ./enh_uni_cluster.bed >./enh_uni_cluseter_SNP_num
#other
  awk '{  if(NR==FNR){ a[$6]++;next }   if(NR>FNR){ if($4 in a){ print $4"\t"a[$4] }else{print $4"\t"0} } }'  ./enh_other_cluster_SNP.bed  ./enh_other_cluster.bed >./enh_other_cluseter_SNP_num
#random
  awk '{  if(NR==FNR){ a[$6]++;next }   if(NR>FNR){ if($4 in a){ print $4"\t"a[$4] }else{print $4"\t"0} } }'  ./randomBed_cluster_SNP.bed  ./randomBed_cluster.bed >./randomBed_cluster_SNP_num

#统计SNP位点距离增强子中心的距离：
cd /media/ding/000B49000006264C/eRNA_project
awk 'BEGIN{OFS="\t"}{ if(NR>1){print $1,$2,$3,$4,$6,$7} }'  /home/ding/all_cmd/script/enh_statistics/enh_all_count_spe_melt | sort -u |sortBed -i - | closest-features --closest  --delim '\t' --dist ./SNP/SNP147_sorted.bed - |awk 'BEGIN{OFS="\t"}($NF<=0&&$NF>=0){print $1,$2,$3,$8,$9}' - >/home/ding/all_cmd/script/enh_statistics/enh_all_SNP.bed

#比较增强子内的SNP距增强子中心的距离：
awk 'BEGIN{OFS="\t"}{ if(NR>1){printf "%s\t%d\t%d\t%s\t%s\t%s\n",$1,($2+$3)/2,($2+$3)/2+1,$4,$6,$7} }'  /home/ding/all_cmd/script/enh_statistics/enh_all_count_spe_melt | sort -u |sortBed -i - | closest-features --closest  --delim '\t' --dist /home/ding/all_cmd/script/enh_statistics/enh_all_SNP.bed - |awk 'BEGIN{OFS="\t"}{ print $4,$5,$10,$11,$12 }' - >/home/ding/all_cmd/script/enh_statistics/enh_all_SNP_dist



#所有的增强子周围的SNP位点：enh_all_count_spe.bed 
cd /home/ding/all_cmd/script/enh_statistics
closest-features --closest  --delim '\t' --dist /media/ding/000B49000006264C/eRNA_project/SNP/SNP147_sorted.bed  ./enh_all_count_spe.bed |awk 'BEGIN{OFS="\t"}{ if($NF==0){print $0} }' - |awk 'BEGIN{OFS="\t";print "SNP_chr","SNP_start","start_end","SNP_name","score","SNP_strand","SNP_observed","SNP_class","SNP_function","enh_chr","enh_start","enh_end","enh_name","is_bid","is_spe","tissue","dist"}{print }' - >/media/ding/000B49000006264C/eRNA_project/figure/table/database/enh_all_SNP
















