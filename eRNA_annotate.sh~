#!/bin/bash

#######################################################################
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
#先进行排序：
sort-bed ./cage_enh/bidir.pairs.bed>./cage_enh/bidir.pairs_sorted.bed
#获取双向转录的转录起始位点：#对双向转录的位点进行过滤，去掉基因内及上游1000bp的：
awk 'BEGIN{FS="\t"}split($11,a,","){print $1"\t"$2+a[1]+1"\t"$3-a[2]"\t"$4}' ./cage_enh/bidir.pairs_sorted.bed |bedops --not-element-of 1 - ../../gencode/protein_coding_gene_up1000 >./cage_enh/bidir.pairs_sorted_TSS.bed

#######################################################################
#各个组织的bid eRNA进行merge,因为cluster的增强子是DHS信息，所以需要合并eRNA:
#临时文件保存到/home/ding/all_cmd/script/enh_statistics/a，因为是不覆盖重定向，所以先要确保不存在这个文件
rm /home/ding/all_cmd/script/enh_statistics/a
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i 
  echo $i"中的双向eRNA聚类"
  awk 'BEGIN{OFS="\t"}{ print $6 }' ./enh_find/enh_bid.bed >>/home/ding/all_cmd/script/enh_statistics/a
done
#整理为bed6格式并排序：
sort -u /home/ding/all_cmd/script/enh_statistics/a | awk 'BEGIN{ FS=":|-";OFS="\t" }{print $1,$2,$3,$0}' - |sortBed -i - > /home/ding/all_cmd/script/enh_statistics/bid_eRNA_cluster.bed

#去掉bid eRNA的中间区域，保留两端转录部分：
bedops -d /home/ding/all_cmd/script/enh_statistics/bid_eRNA_cluster.bed   /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/cage_enh/bidir.pairs_sorted_TSS.bed >/home/ding/all_cmd/script/enh_statistics/a

#添加转录方向，左边的是-，右边的是+,并控制范围是TSS上下游500bp: $3-400,$3  $2,$2+400
awk 'BEGIN{OFS="\t"}{ if(NR%2==1){print $1,$3-400,$3,"bid_eRNA_"NR,".\t-"} if(NR%2==0){print $1,$2,$2+400,"bid_eRNA_"NR,".\t+"} }' /home/ding/all_cmd/script/enh_statistics/a > /home/ding/all_cmd/script/enh_statistics/bid_eRNA_cluster_transcript.bed

rm /home/ding/all_cmd/script/enh_statistics/a

#######################################################################
#各个组织的unbid eRNA进行merge
#临时文件保存到/home/ding/all_cmd/script/enh_statistics/a，因为是不覆盖重定向，所以先要确保不存在这个文件
rm /home/ding/all_cmd/script/enh_statistics/a
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i 
  echo $i"中的单向eRNA聚类"
   awk 'BEGIN{OFS="\t"}{ print $6 }' ./enh_find_1/enh_unbid.bed >>/home/ding/all_cmd/script/enh_statistics/a
done
#整理为bed6格式并排序然后取TSS后400bp：$3-400,$3  $2,$2+400
sort -u /home/ding/all_cmd/script/enh_statistics/a |awk 'BEGIN{FS=":|-|,";OFS="\t"}{ print $1,$2,$3,"unbid_eRNA_"NR,".",$4 }' -  |awk 'BEGIN{OFS="\t"}{ if($6!="+"){print $1,$2,$3,$4,$5,"-"} if($6=="+")print $0 }' - |sortBed -i - |awk 'BEGIN{FS="\t";OFS="\t"}{if($NF=="-"){print $1,$3-400,$3,$4,$5,$NF} if($NF=="+"){print $1,$2,$2+400,$4,$5,$NF}}' - > /home/ding/all_cmd/script/enh_statistics/unbid_eRNA_cluster_transcript.bed

rm /home/ding/all_cmd/script/enh_statistics/a

#合并unbid和bid的eRNA到一个文件：
cd /home/ding/all_cmd/script/enh_statistics
awk '{print }' ./bid_eRNA_cluster_transcript.bed ./unbid_eRNA_cluster_transcript.bed >./all_eRNA_cluster_transcript.bed

#######################################################################
#与除去mRNA的注释文件TSS比对，与TSS+400有交集即算：/media/ding/000B49000006264C/eRNA_project/gencode/TSS_no_mRNA.bed
cd /home/ding/all_cmd/script/enh_statistics
closestBed -s -d  -a ./all_eRNA_cluster_transcript.bed -b /media/ding/000B49000006264C/eRNA_project/gencode/TSS_no_mRNA.bed | awk '{ if($NF==0)print }' - >./eRNA_annotate













