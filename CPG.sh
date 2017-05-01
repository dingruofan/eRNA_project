#!/bin/bash

#获得所有enh文件的注释：
#增强子中心与CpG岛的距离，
awk 'BEGIN{FS="\t";OFS="\t";}{ if(NR>1){print $1,($2+$3)/2,($2+$3)/2+1,$4,$5,$6} }' /home/ding/all_cmd/script/enh_statistics/enh_all|sortBed -i - |closest-features --closest  --delim '\t' --dist  - /media/ding/000B49000006264C/eRNA_project/CPG/CPG.bed|awk 'BEGIN{FS="\t";OFS="\t"}{ if($NF!="NA"){print $4,$5,$6,$NF}  }' - >/home/ding/all_cmd/script/enh_statistics/CPG_dist



closest-features --closest  --delim '\t' --dist /media/ding/000B49000006264C/eRNA_project/gencode/protein_gene2.ss /media/ding/000B49000006264C/eRNA_project/CPG/CPG.bed|awk 'BEGIN{FS="\t";OFS="\t"}{ if($NF!="NA"){print "gene_"NR,"type","isspe",$NF}  }'>/home/ding/all_cmd/script/enh_statistics/CPG_protein_gene2_dist


#CPG与基因TSS的关系：
closest-features --closest  --delim '\t' --dist /media/ding/000B49000006264C/eRNA_project/CPG/CPG.bed  /media/ding/000B49000006264C/eRNA_project/gencode/protein_gene2.ss|awk 'BEGIN{FS="\t";OFS="\t"}{ if($NF!="NA"&&NF>-1000&&NF<1000){print }  }' |wc -l 


#CPG与所有gencode gene TSS的关系：
closest-features --closest  --delim '\t' --dist /media/ding/000B49000006264C/eRNA_project/CPG/CPG.bed  /media/ding/000B49000006264C/eRNA_project/gencode/gencode.v19.annotation.TSS |awk 'BEGIN{FS="\t";OFS="\t"}{ if($NF!="NA"&&NF>-1000&&NF<1000){print }  }' |wc -l 


closest-features --closest  --delim '\t' --dist /media/ding/000B49000006264C/eRNA_project/CPG/CPG.bed  /media/ding/000B49000006264C/eRNA_project/gencode/gencode.v19.annotation.TSS |awk 'BEGIN{FS="\t";OFS="\t"}{ if($NF!="NA"&&NF>-1000&&NF<1000){print }  }' |wc -l 



##################################
#增强子与CpG岛的距离：
################################
#获得所有enh文件的注释：
awk 'BEGIN{FS="\t";OFS="\t";}{ if(NR>1){print $1,$2,$3,$4,$5,$6} }' /home/ding/all_cmd/script/enh_statistics/enh_all|sortBed -i - |closest-features --closest  --delim '\t' --dist  - /media/ding/000B49000006264C/eRNA_project/CPG/CPG.bed|awk 'BEGIN{FS="\t";OFS="\t"}{ if($NF!="NA"){print $4,$5,$6,$NF}  }' - >/home/ding/all_cmd/script/enh_statistics/CPG_dist





closest-features --closest  --delim '\t' --dist /media/ding/000B49000006264C/eRNA_project/gencode/protein_gene2.ss /media/ding/000B49000006264C/eRNA_project/CPG/CPG.bed|awk 'BEGIN{FS="\t";OFS="\t"}{ if($NF!="NA"){print "gene_"NR,"type","isspe",$NF}  }'>/home/ding/all_cmd/script/enh_statistics/CPG_protein_gene2_dist








