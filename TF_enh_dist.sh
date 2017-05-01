#!/bin/bash

#######################################################################
#操作enh_all_count文件和/media/ding/000B49000006264C/eRNA_project/TFBS/TFBS_sorted.bed
cd /home/ding/all_cmd/script/enh_statistics
#这个是按照bid/unbid/no_eRNA的enh进行分类的：
awk 'BEGIN{OFS="\t"}{ if(NR>1){printf "%s\t%d\t%d\t%s\n",$1,($2+$3)/2,($2+$3)/2+1,$NF} }' enh_all_count |sortBed -i -| closest-features --closest  --delim '\t' --dist  /media/ding/000B49000006264C/eRNA_project/TFBS/TFBS_sorted.bed - |awk 'BEGIN{OFS="\t"}{ if($NF>-10000&&$NF<10000){print $4,$9,$NF} }' - >./TF_enh_dist_iseRNA_type

#这个是那照spe/uni/other的enh进行分类的：
cd /home/ding/all_cmd/script/enh_statistics
awk 'BEGIN{OFS="\t"}{ if(NR>1){printf "%s\t%d\t%d\t%s\n",$1,($2+$3)/2,($2+$3)/2+1,$NF} }' enh_spe |sortBed -i -|closest-features --closest  --delim '\t' --dist  /media/ding/000B49000006264C/eRNA_project/TFBS/TFBS_sorted.bed - |awk 'BEGIN{OFS="\t"}{ if($NF>-10000&&$NF<10000){print $4,$9,$NF} }' ->./TF_enh_dist_isspe_type

echo "数据准备完成"
#使用R进行分析：













