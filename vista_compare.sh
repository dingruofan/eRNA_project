#!/bin/bash

cd /home/ding/all_cmd/script
#########################################

#统计vista enhancer的长度，选择基因外的vista enhancer进行比较：
bedops --not-element-of 1 /home/ding/all_cmd/script/enh_statistics/vista_enhancer.bed /media/ding/000B49000006264C/eRNA_project/gencode/protein_coding_gene_up1000 > /home/ding/all_cmd/script/enh_statistics/vista_enhancer_outgene.bed

#vista enh平均长度达到了1783,所以先对我得到的enh进行扩增长度：以enh中心上下游各扩增1000bp:
awk 'BEGIN{OFS="\t"}{ if(NR>1){printf "%s\t%d\t%d\t%s\n",$1,$2,$3,$NF} }' /home/ding/all_cmd/script/enh_statistics/enh_all_count |sortBed -i -| closest-features --closest  --delim '\t' --dist /home/ding/all_cmd/script/enh_statistics/vista_enhancer_outgene.bed  - |awk '{print $7,$NF}'  - >/home/ding/all_cmd/script/enh_statistics/dist_to_plot


#与vista enh脑的数据比较：forebrain
bedops --not-element-of 1 /media/ding/000B49000006264C/eRNA_project/vista_enhancer/forebrain/vista_enh.bed /media/ding/000B49000006264C/eRNA_project/gencode/protein_coding_gene_up1000 |sort-bed - > /media/ding/000B49000006264C/eRNA_project/vista_enhancer/forebrain/vista_enh_outgene.bed
#提取brain的enh,第六列：
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$6==1){printf "%s\t%d\t%d\t%s\n",$1,$2,$3,$NF} }' /home/ding/all_cmd/script/enh_statistics/enh_all_count |sortBed -i -| closest-features --closest  --delim '\t' --dist /media/ding/000B49000006264C/eRNA_project/vista_enhancer/forebrain/vista_enh_outgene.bed  - |awk '{print $7,$NF}'  - >/home/ding/all_cmd/script/enh_statistics/vista_brain_dist_to_plot

#与vista enh脑的数据比较：midbrain
bedops --not-element-of 1 /media/ding/000B49000006264C/eRNA_project/vista_enhancer/midbrain_mesencephalon/vista_enh.bed /media/ding/000B49000006264C/eRNA_project/gencode/protein_coding_gene_up1000 |sort-bed - > /media/ding/000B49000006264C/eRNA_project/vista_enhancer/midbrain_mesencephalon/vista_enh_outgene.bed
#提取brain的enh,第六列：
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$6==1){printf "%s\t%d\t%d\t%s\n",$1,$2,$3,$NF} }' /home/ding/all_cmd/script/enh_statistics/enh_all_count |sortBed -i -| closest-features --closest  --delim '\t' --dist /media/ding/000B49000006264C/eRNA_project/vista_enhancer/midbrain_mesencephalon/vista_enh_outgene.bed  - |awk '{print $7,$NF}'  - >/home/ding/all_cmd/script/enh_statistics/vista_brain_dist_to_plot

####################这个也要,包括midbrain和forebrain：
#基因外的：
#提取brain的enh,第六列：
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$6==1){printf "%s\t%d\t%d\t%s\n",$1,$2,$3,$NF} }' /home/ding/all_cmd/script/enh_statistics/enh_all_count |sortBed -i -| closest-features --closest  --delim '\t' --dist /media/ding/000B49000006264C/eRNA_project/vista_enhancer/fore_mid_brain/vista_enh_outgene.bed  - |awk '{print $7,$NF}'  - >/home/ding/all_cmd/script/enh_statistics/vista_brain_dist_to_plot

##################


#############这个是要的
#与vista enh心脏的数据比较：
bedops --not-element-of 1 /media/ding/000B49000006264C/eRNA_project/vista_enhancer/heart/vista_enh.bed /media/ding/000B49000006264C/eRNA_project/gencode/protein_coding_gene_up1000 > /media/ding/000B49000006264C/eRNA_project/vista_enhancer/heart/vista_enh_outgene.bed
#提取heart的enh,第八列：
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$8==1){printf "%s\t%d\t%d\t%s\n",$1,$2,$3,$NF} }' /home/ding/all_cmd/script/enh_statistics/enh_all_count |sortBed -i -| closest-features --closest  --delim '\t' --dist /media/ding/000B49000006264C/eRNA_project/vista_enhancer/heart/vista_enh_outgene.bed  - |awk '{print $7,$NF}'  - >/home/ding/all_cmd/script/enh_statistics/vista_heart_dist_to_plot
#############


#与heart一致的：
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$8==1){printf "%s\t%d\t%d\t%s\n",$1,$2,$3,$4,$NF} }' /home/ding/all_cmd/script/enh_statistics/enh_all_count |sortBed -i -| closest-features --closest  --delim '\t' --dist /media/ding/000B49000006264C/eRNA_project/vista_enhancer/heart/vista_enh_outgene.bed  - |awk 'BEGIN{OFS="\t"}{if($NF==0)print $4,$5,$6,$7}' - > ~/桌面/consistent_heart_vista.bed

#与brain一致的：
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$8==1){printf "%s\t%d\t%d\t%s\n",$1,$2,$3,$4,$NF} }' /home/ding/all_cmd/script/enh_statistics/enh_all_count |sortBed -i -| closest-features --closest  --delim '\t' --dist /media/ding/000B49000006264C/eRNA_project/vista_enhancer/fore_mid_brain/vista_enh_outgene.bed  - |awk 'BEGIN{OFS="\t"}{if($NF==0)print $4,$5,$6,$7}' - > ~/桌面/consistent_brain_vista.bed




#与vista enh liver的数据比较：
bedops --not-element-of 1 /media/ding/000B49000006264C/eRNA_project/vista_enhancer/liver/vista_enh.bed /media/ding/000B49000006264C/eRNA_project/gencode/protein_coding_gene_up1000 > /media/ding/000B49000006264C/eRNA_project/vista_enhancer/liver/vista_enh_outgene.bed
#提取brain的enh,第九列：
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$9==1){printf "%s\t%d\t%d\t%s\n",$1,($2+$3)/2-1000,($2+$3)/2+1000,$NF} }' /home/ding/all_cmd/script/enh_statistics/enh_all_count |sortBed -i -| closest-features --closest  --delim '\t' --dist /media/ding/000B49000006264C/eRNA_project/vista_enhancer/liver/vista_enh_outgene.bed  - |awk '{print $7,$NF}'  - >/home/ding/all_cmd/script/enh_statistics/vista_liver_dist_to_plot


cd /home/ding/all_cmd/script
now_path=`pwd`
Rscript /home/ding/all_cmd/script/vista_compare.R  $now_path
feh $now_path/result.mean_signal.png





#做vista enhancer中心周围的GC含量图：
cd /home/ding/all_cmd/script/
bwtool agg 2000:2000 ./enh_statistics/vista_enhancer_outgene.bed,./enh_statistics/random_all_enh.bed,./enh_statistics/enh_bid_cluster.bed,./enh_statistics/enh_unbid_cluster.bed,/media/ding/000B49000006264C/eRNA_project/gencode/protein_gene2.ss,./enh_statistics/no_eRNA_enh.bed /media/ding/000B49000006264C/eRNA_project/GC/hg19.gc5Base.bw  /dev/stdout > ./enh_statistics/GC_vista_result.mean_signal









