#!/bin/bash

#######################################################################
#新思路2：筛选Dnase位点2000bp范围内组蛋白修饰的峰值summits，然而并没有summits
#######################################################################

#2kb内有me1/27ac:
closest-features --closest  --delim '\t' --dist ../../DNase/$1/Dnase_ENCODE_cluster.bed ./bed/H3K27ac/H3K27ac_peaks.bed|awk '($NF<=2000&&$NF>=-2000){print $1,$2,$3}' - |sed  's/ /\t/g' - >./enh_find/a

closest-features --closest  --delim '\t' --dist ./enh_find/a ./bed/H3K4me1/H3K4me1_peaks.bed |awk '($NF<=2000&&$NF>=-2000){print $1,$2,$3}' - |sed  's/ /\t/g' - >./enh_find/b


#去掉蛋白编码区：
bedops --not-element-of 1 ./enh_find/b ../../gencode/protein_coding_gene_up1000 >./enh_find/c

#H3K4me3的处理：然而并没有处理：($NF>0||$NF<-0)($NF>=2000||$NF<=-2000)
closest-features --closest  --delim '\t' --dist ./enh_find/c ./bed/H3K4me3/H3K4me3_peaks.bed |awk '{print $1,$2,$3}' - |sed  's/ /\t/g' - >./enh_find/result.bed


#处理信号值：要求enh中心附近2kb内H3K4me1的信号的均值大于H3K4me3的信号的均值，没有信号的位置不做统计。(第八列是均值) -header参数查看列名。
bwtool summary ./enh_find/result.bed ./bg/H3K4me1.bw -with-sum  /dev/stdout |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$8}' ->./enh_find/H3K4me1_sig_2kb_mean
bwtool summary ./enh_find/result.bed ./bg/H3K4me3.bw -with-sum  /dev/stdout |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$8}' ->./enh_find/H3K4me3_sig_2kb_mean
#进行比较：
awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1"\t"$2"\t"$3]=$4}  if(NR>FNR&&a[$1"\t"$2"\t"$3]-$4>0){print $1,$2,$3} }' ./enh_find/H3K4me1_sig_2kb_mean ./enh_find/H3K4me3_sig_2kb_mean |sortBed -i - >./enh_find/result.bed

rm ./enh_find/H3K4me1_sig_2kb_mean ./enh_find/H3K4me3_sig_2kb_mean

#做蛋白信号图：
bwtool agg 2000:2000 ./enh_find/result.bed ./bg/H3K4me1.bw,./bg/H3K4me3.bw,./bg/H3K27ac.bw,../../DNase/$1/new_DHS.bigWig  /dev/stdout > ./enh_find/result.mean_signal
now_path=`pwd`
Rscript $2/bid_plot.R $now_path
echo $1
feh $now_path/result.mean_signal.png
rm $now_path/result.mean_signal.png


#######################################################################
#与CAGE定义的双向转录取交集：
#组蛋白定义的enhancer中心result.bed上下游扩增2kb:

awk '{print $1"\t"$2-2000"\t"$3+2000}' ./enh_find/result.bed  >./enh_find/result_2kb.bed

#2kb内有双向转录TSS的增强子：(bidir.pairs_sorted_TSS.bed)
#一个dnase 2kb范围内可能存在多个bidir，选择离他最近的那一个。
bedops --element-of 100% ../../fantom/$1_CAGE/cage_enh/bidir.pairs_sorted_TSS.bed ./enh_find/result_2kb.bed  |closest-features --closest  --delim '\t' --dist -  ./enh_find/result.bed |awk '{print $5"\t"$6"\t"$7"\t"$4"\t.\t"$NF}' - >./enh_find/result_2kb_bidir_all_1.bed 


#下面这条命令是：DHS 2000bp范围内只要有双向转录的既为合格的
bedops --element-of 100% ../../fantom/$1_CAGE/cage_enh/bidir.pairs_sorted_TSS.bed ./enh_find/result_2kb.bed  |closest-features --closest  --delim '\t' --dist  ./enh_find/result.bed - | awk '($NF<=2000&&$NF>=-2000){print $0}' - | awk '{print $1"\t"$2"\t"$3"\t"$7"\t.\t"$NF}' - >./enh_find/result_2kb_bidir_all_2.bed 
 

#上面两条取交集就可以得到无重复的： 
awk '{ if(NR==FNR){a[$1"\t"$2"\t"$3"\t"$4]} if(NR>FNR&&($1"\t"$2"\t"$3"\t"$4 in a)){print $0} }' ./enh_find/result_2kb_bidir_all_2.bed ./enh_find/result_2kb_bidir_all_1.bed >./enh_find/result_2kb_bidir_all.bed 

#进行排序：
sort-bed ./enh_find/result_2kb_bidir_all.bed >./enh_find/Dnase_histone_bidir_enh.bed

#到此得到的./enh_find/result_2kb_bidir_all.bed 包含cage定义的双向转录为点（前4列，第四列是完整的），以及5,6,7列是与之最近的2000bp内的Dnase位点。

#处理一下列的位置,并注意需要重新排序：
#awk '{print $5"\t"$6"\t"$7"\t"$4"\t.\t"$NF}'  ./enh_find/result_2kb_bidir_all.bed |sort-bed - > ./enh_find/Dnase_histone_bidir_enh.bed
#awk '{print $1"\t"$2"\t"$3"\t"$7"\t.\t"$NF}'  ./enh_find/result_2kb_bidir_all.bed |sort-bed - > ./enh_find/Dnase_histone_bidir_enh.bed


#######################################################################

#获得2000bp内的转录因子及其结合位点。(注意有些enh附近没有一个启动子)

closest-features --closest  --delim '\t' --dist ../../TFBS/TFBS_sorted.bed ./enh_find/Dnase_histone_bidir_enh.bed |awk '($NF<=2000&&$NF>=-2000){print }' - |sed  's/ /\t/g' - >./enh_find/Dnase_histone_bidir_enh_TFBS

#列说明：前5列为转录因子结合为点信息,6,7,8列为Dnase结合为点，9列为bidir位点，10列为strand,11列为Dnase与bidir距离，12列为Dnase与TFBS距离。

#######################################################################
#添加TF表大量到最后一列：(表达量文件../../TF_exp/$1/TF_exp)
awk 'BEGIN{FS="\t"}{if(NR==FNR){a[$2]=$3;next} if(NR>FNR){print $0"\t"a[$4]}}'  ../../TF_exp/$1/TF_exp ./enh_find/Dnase_histone_bidir_enh_TFBS >./enh_find/Dnase_histone_bidir_enh_TFBS_exp 


#处理一下列的位置，并注意需要排序：
awk '{print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$1"\t"$2"\t"$3"\t"$4"\t"$NF}' ./enh_find/Dnase_histone_bidir_enh_TFBS_exp |sort-bed -> ./enh_find/Dnase_histone_bidir_enh_TFBS_exp.bed


#将统计信息写入文件：
cat >>$2/enh_statistics/bid_statistics<<EOF
$1	`awk 'END{print NR}' ./enh_find/result.bed`	`awk 'END{print NR}' ./enh_find/result_2kb_bidir_all.bed`	`awk 'END{print NR}' ./enh_find/Dnase_histone_bidir_enh_TFBS_exp.bed` 
EOF

echo -n "符合组蛋白修饰的有："
wc -l ./enh_find/result.bed
echo -n "enh共有："
wc -l ./enh_find/Dnase_histone_bidir_enh.bed 
echo -n "转录因子为点有："
wc -l ./enh_find/Dnase_histone_bidir_enh_TFBS_exp.bed


