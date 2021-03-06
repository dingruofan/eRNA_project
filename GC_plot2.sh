#!/bin/bash
#使用的是没有strand的TSS：
cd /home/ding/all_cmd/script/
n=`wc -l ./enh_statistics/enh_all_count|awk '{print $1-1}' - `
randomBed -l 200 -n $n -seed 2 -g ../hg19.chrom_24.sizes|sortBed -i - |bedops --not-element-of 1 - /media/ding/000B49000006264C/eRNA_project/gencode/protein_coding_gene.bed >./enh_statistics/random_all_enh.bed

bwtool agg 2000:2000 ./enh_statistics/random_all_enh.bed,./enh_statistics/enh_spe_cluster.bed,./enh_statistics/enh_uni_cluster.bed,/media/ding/000B49000006264C/eRNA_project/gencode/protein_gene2.ss,./enh_statistics/enh_other_cluster.bed /media/ding/000B49000006264C/eRNA_project/GC/hg19.gc5Base.bw  /dev/stdout > ./enh_statistics/GC2_result.mean_signal



