#!/bin/bash
#######################################################################
#处理10个组织的CAGE，即合并TCs后进行分析双向转录：
mkdir /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
/home/ding/tool/enhancers-master/scripts/bidir_enhancers -d 0.8  -f ../ctss.path -o ./cage_enh  

#由于CAGE_master生成的是双向转录的TPM，而我们需要bid转录的每个的表达量，所以先用bedtools的bed12ToBed6将双向的转为bed6，这样就分割了双向到两行，然后再转为bed12格式。
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
bed12ToBed6 -i ./cage_enh/bidir.pairs.bed | awk 'BEGIN{OFS="\t";} { if(NR%2==1){print $1,$2,$3,$4",-",$5,"-"} if(NR%2==0){print $1,$2,$3,$4",+",$5,"+"} } ' - |sortBed -i - >./cage_enh/bidir_pairs_sorted.bed6
#转为bed12格式：
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
awk '{ if($1=="chrM"){next}else{ print $0"\t"$2"\t"$3"\t0,0,0\t1\t"$3-$2"\t0"} }' ./cage_enh/bidir_pairs_sorted.bed6 >./cage_enh/bidir_pairs_sorted.bed12
#使用quantify_enhancers程序获得每个双向转录本的的表达量。
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
/home/ding/tool/enhancers-master/scripts/quantify_enhancers -f ../ctss.path -e ./cage_enh/bidir_pairs_sorted.bed12 -o ./bid_cage_enh_quantify
#因为程序有bug,导致有空行，但不影响，提取表达量文件并进行归一化。
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
awk '{if(NF>5)print $0}' ./bid_cage_enh_quantify/enhancers.expression.matrix >./bid_cage_enh_quantify/bidir_pairs_sorted_exp.matrix
awk '{if(NF>5)print $0}' ./bid_cage_enh_quantify/enhancers.expression.tpm.matrix >./bid_cage_enh_quantify/bidir_pairs_sorted_exp_tpm.matrix
#对TPM.Matrix文件添加列名：
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
awk 'BEGIN{print "genome_location\tAdrenal\tBrain\tBreast\tHeart\tLiver\tLung\tOvary\tPlacenta\tSkeletalMuscle\tKidney\tAdipose\tPancreas"}{print $0}' ./cage_enh/bidir.pairs.expression.tpm.matrix>./cage_enh/bidir_pairs_expression_tpm.matrix
awk 'BEGIN{print "genome_location\tAdrenal\tBrain\tBreast\tHeart\tLiver\tLung\tOvary\tPlacenta\tSkeletalMuscle\tKidney\tAdipose\tPancreas"}{print $0}' ./bid_cage_enh_quantify/bidir_pairs_sorted_exp_tpm.matrix >./bid_cage_enh_quantify/bidir_pairs_expression_tpm.matrix

#注意：bid_cage_enh_quantify目录下放的marix是后来生成的表达量文件，而cage_enh下放的matrix是最初统一生成的matrix文件。
#在统计enh的时候用的是cage_enh下的，在计算表达量的时候用的是后面生成的。




###################

#对TCs进行排序并生成bed12格式的TCs，只有bed12才可以使用quantify_enhancers程序
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
sort-bed ./cage_enh/TCs.bed>./cage_enh/TCs_sorted.bed
awk '{ if($1=="chrM"){next}else{ print $0"\t"$2"\t"$3"\t0,0,0\t1\t"$3-$2"\t0"} }' ./cage_enh/TCs_sorted.bed>./cage_enh/TCs_sorted.bed12

#使用quantify_enhancers程序获得每个TCs的表达量。
/home/ding/tool/enhancers-master/scripts/quantify_enhancers -f ../ctss.path -e ./cage_enh/TCs_sorted.bed12 -o ./cage_enh_quantify

#因为程序有bug,导致有空行，但不影响，提取表达量文件并进行归一化。
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
awk '{if(NF>5)print $0}' ./cage_enh_quantify/enhancers.expression.matrix >./cage_enh_quantify/TCs_sorted_exp.matrix
awk '{if(NF>5)print $0}' ./cage_enh_quantify/enhancers.expression.tpm.matrix >./cage_enh_quantify/TCs_sorted_exp_tpm.matrix
#rm  ./cage_enh_quantify/* ./cage_enh/TCs_sorted.bed12 ./cage_enh/TCs.bed12
#对TPM.Matrix文件添加列名：
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
awk 'BEGIN{print "genome_location\tAdrenal\tBrain\tBreast\tHeart\tLiver\tLung\tOvary\tPlacenta\tSkeletalMuscle\tKidney\tAdipose\tPancreas"}{print $0}'  ./cage_enh_quantify/TCs_sorted_exp_tpm.matrix >./cage_enh/TCs_sorted_exp_tpm.matrix



###############################################################################
#处理双端转录的Matrix，并将bidir.pairs.bed提取到每一个组织中，表达量所对列分别是：2-11
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney Adipose Pancreas)
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
for i in ${tissue[@]}
do
  awk -v tissue="$i" 'BEGIN{ j=0 } { if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){j=i}} next} if(NR==FNR&&NR>1&&$j!=0){a[$1]=$j} if(NR>FNR&&($4 in a)){print $1"\t"$2"\t"$3"\t"$4"\t"a[$4]"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12} }' ./cage_enh/bidir_pairs_expression_tpm.matrix ./cage_enh/bidir.pairs.bed  > ../$i"_CAGE"/cage_enh/bidir.pairs.bed
done


#处理单端转录的Matrix，并将TCs.bed提取到每一个组织中：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney  Adipose Pancreas)
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
for i in ${tissue[@]}
do
  awk -v tissue="$i" 'BEGIN{ j=0 } { if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){j=i}} next} if(NR==FNR&&NR>1&&$j!=0){a[$1]=$j} if(NR>FNR&&($4 in a)){print $1"\t"$2"\t"$3"\t"$4"\t"a[$4]"\t"$6} }' ./cage_enh/TCs_sorted_exp_tpm.matrix ./cage_enh/TCs_sorted.bed > ../$i"_CAGE"/cage_enh/TCs.bed
done 


#进行排序和去掉基因内的TCs和bidir.pairs。
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney Adipose Pancreas)
for i in ${tissue[@]}
do
  cd /media/ding/000B49000006264C/eRNA_project/fantom/$i"_CAGE"
  #先进行排序：
  sort-bed ./cage_enh/enhancers.bed>./cage_enh/enhancers_sorted.bed
  sort-bed ./cage_enh/bidir.pairs.bed>./cage_enh/bidir.pairs_sorted.bed
  sort-bed ./cage_enh/TCs.bed> ./cage_enh/TCs_sorted.bed
  
  
  #获取双向转录的转录起始位点：#对双向转录的位点进行过滤，去掉基因内及上游1000bp的：
  awk 'BEGIN{FS="\t"}split($11,a,","){print $1"\t"$2+a[1]"\t"$3-a[2]"\t"$4}' ./cage_enh/bidir.pairs_sorted.bed |bedops --not-element-of 1 - ../../gencode/protein_coding_gene_up1000 >./cage_enh/bidir.pairs_sorted_TSS.bed

  #去掉双端转录和基因以及基因上游1000bp内包含的TCs，获得基因外单向转录的TCs
  bedops --not-element-of 1 ./cage_enh/TCs_sorted.bed ./cage_enh/bidir.pairs_sorted.bed |bedops --not-element-of 1 - ../../gencode/protein_coding_gene_up1000 >./cage_enh/TCs_sorted_extragene.bed

done





############################################################
#在这里将表达量TPM转为BPKM：
#双向的：TPM/10^6*all_count/transcript_length
#提取双向转录本两个转录本长度的和：
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/cage_enh
awk 'BEGIN{OFS="\t";} { split($11,a,","); print $1,$2,$3,$4,a[1]+a[2]} ' bidir.pairs.bed >bidir.pairs_transc_length.bed

#提取单向转录本的长度：
awk 'BEGIN{OFS="\t";} { print $1,$2,$3,$4,$11 } ' TCs_sorted.bed12 >TCs_sorted_transc_length.bed

#接下来使用R脚本比较方便：
cd /home/ding/all_cmd/script
Rscript CAGE_all_tissue.R

#













