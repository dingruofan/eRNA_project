#!/bin/bash

#######################################################################
#基因转录起始位点(/media/ding/000B49000006264C/eRNA_project/gencode/all_gene_protein_gene_TSS_sorted.bed)
#需要对正负链进行调整：
cd /home/ding/all_cmd/script/enh_statistics
closest-features --closest  --delim '\t' --dist ./enh_bid_cluster.bed /media/ding/000B49000006264C/eRNA_project/gencode/all_gene_protein_gene_TSS_sorted.bed | awk 'BEGIN{OFS="\t"} { if($2>$6&&$10=="+"&&$NF<0){print $0,-1*$NF;}else if($2<$6&&$10=="+"&&$NF>0){print $0,-1*$NF;}else if($2>$6&&$10!="-"&&$NF>0){print $0,-1*$NF;} else if($2<$6&&$10!="+"&&$NF<0){print $0,-1*$NF;}else{print $0,$NF} }' - >./a

closest-features --closest  --delim '\t' --dist ./enh_unbid_cluster.bed /media/ding/000B49000006264C/eRNA_project/gencode/all_gene_protein_gene_TSS_sorted.bed | awk 'BEGIN{OFS="\t"} { if($2>$6&&$10=="+"&&$NF<0){print $0,-1*$NF;}else if($2<$6&&$10=="+"&&$NF>0){print $0,-1*$NF;}else if($2>$6&&$10!="-"&&$NF>0){print $0,-1*$NF;} else if($2<$6&&$10!="+"&&$NF<0){print $0,-1*$NF;}else{print $0,$NF} }' - >./b

closest-features --closest  --delim '\t' --dist ./enh_no_eRNA_cluster.bed /media/ding/000B49000006264C/eRNA_project/gencode/all_gene_protein_gene_TSS_sorted.bed | awk 'BEGIN{OFS="\t"} { if($2>$6&&$10=="+"&&$NF<0){print $0,-1*$NF;}else if($2<$6&&$10=="+"&&$NF>0){print $0,-1*$NF;}else if($2>$6&&$10!="-"&&$NF>0){print $0,-1*$NF;} else if($2<$6&&$10!="+"&&$NF<0){print $0,-1*$NF;}else{print $0,$NF} }' - >./c

awk 'BEGIN{OFS="\t"}{ if(ARGIND==1){print "bid",$NF} if(ARGIND==2){print "unbid",$NF} if(ARGIND==3){print "no_eRNA",$NF} }' ./a ./b ./c >./enh_protein_TSS_dist

rm ./a ./b ./c


#spe的：
cd /home/ding/all_cmd/script/enh_statistics
closest-features --closest  --delim '\t' --dist ./enh_spe_cluster.bed /media/ding/000B49000006264C/eRNA_project/gencode/all_gene_protein_gene_TSS_sorted.bed | awk 'BEGIN{OFS="\t"} { if($2>$6&&$10=="+"&&$NF<0){print $0,-1*$NF;}else if($2<$6&&$10=="+"&&$NF>0){print $0,-1*$NF;}else if($2>$6&&$10!="-"&&$NF>0){print $0,-1*$NF;} else if($2<$6&&$10!="+"&&$NF<0){print $0,-1*$NF;}else{print $0,$NF} }' - >./a

closest-features --closest  --delim '\t' --dist ./enh_other_cluster.bed /media/ding/000B49000006264C/eRNA_project/gencode/all_gene_protein_gene_TSS_sorted.bed | awk 'BEGIN{OFS="\t"} { if($2>$6&&$10=="+"&&$NF<0){print $0,-1*$NF;}else if($2<$6&&$10=="+"&&$NF>0){print $0,-1*$NF;}else if($2>$6&&$10!="-"&&$NF>0){print $0,-1*$NF;} else if($2<$6&&$10!="+"&&$NF<0){print $0,-1*$NF;}else{print $0,$NF} }' - >./b

closest-features --closest  --delim '\t' --dist ./enh_uni_cluster.bed /media/ding/000B49000006264C/eRNA_project/gencode/all_gene_protein_gene_TSS_sorted.bed | awk 'BEGIN{OFS="\t"} { if($2>$6&&$10=="+"&&$NF<0){print $0,-1*$NF;}else if($2<$6&&$10=="+"&&$NF>0){print $0,-1*$NF;}else if($2>$6&&$10!="-"&&$NF>0){print $0,-1*$NF;} else if($2<$6&&$10!="+"&&$NF<0){print $0,-1*$NF;}else{print $0,$NF} }' - >./c

awk 'BEGIN{OFS="\t"}{ if(ARGIND==1){print "spe",$NF} if(ARGIND==2){print "oth",$NF} if(ARGIND==3){print "uni",$NF} }' ./a ./b ./c >./enh_protein_TSS_dist2

rm ./a ./b ./c










