#!/bin/bash

#######################################################################
#处理从tiger得到的特异性gene，合并成一个文件：
#注意：mammary_gland用来代替breast，因为没有breast。muscle代替skeletamuscle：

#####################################
#下面使用pagenbase的：

#每个组织的TS enh/TF/target 的数目分布:
# 
# cd /media/ding/000B49000006264C/eRNA_project/tiger
# #创建结果文件：
# :>/home/ding/all_cmd/script/enh_statistics/TS_gene
# awk 'BEGIN{print "gene_name\ttissue"}' ./TS_gene > /home/ding/all_cmd/script/enh_statistics/TS_gene
# tissue=(Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney )
# #tissue=(Adrenal)
# for i in ${tissue[@]}
# do
# awk -v tissue=$i 'BEGIN{OFS="\t"}{ if(NR>2){print $3,tissue} }' $i"_genes.txt"|sort -u >>/home/ding/all_cmd/script/enh_statistics/TS_gene
# echo $i
# done 

#与enh_exp_targetgene_exp（由./enh_target_gene_exp.sh产生）比较，获得enh调控的特异性的靶基因，(并且靶基因的表达量>0:先不用)
cd /home/ding/all_cmd/script/enh_statistics/
#awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1"\t"$3]=$2} if((NR>FNR)&&($10"\t"$13 in a)&&($12!="random")&&($11>0)){print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,a[$10"\t"$13]} }' ./pagenbase_ts_gene ./enh_exp_targetgene_exp >./enh_exp_TStargetgene_exp
awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1"\t"$3]=$2} if((NR>FNR)&&($10"\t"$13 in a)&&($12!="random")){print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,a[$10"\t"$13]} }' ./pagenbase_ts_gene ./enh_exp_targetgene_exp >./enh_exp_TStargetgene_exp

#生成的enh_exp_TStargetgene_exp与spe_TF_enh_tissue_3（由enh_ts.sh生成）比较,得到TS enh对应的TS TF和对应的TS gene:
cd /home/ding/all_cmd/script/enh_statistics/
awk 'BEGIN{OFS="\t"}{ print $4,$12,$9,$8  }' ./enh_exp_TStargetgene_exp |sort -u | awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1"\t"$2]=$3"\t"$4} if(NR>FNR&&($1"\t"$4 in a)){print $0,a[$1"\t"$4]} }' - ./spe_TF_enh_tissue_3 >./TS_enh_TF_targetgene


#特异性enh与特异性的targetgene,(并且target gene的表达量>0,先不用)
#处理上面得到的：./enh_exp_TStargetgene_exp与./enh_spe，结果文件中两列数字分别是gene的表达量和特异性分数。
cd /home/ding/all_cmd/script/enh_statistics/
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$NF=="spe"){print $4,$NF} }' ./enh_spe | awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1]=$2} if(NR>FNR&&$4 in a&&$10>0){print $4,$8,$9,$10,$11,$12,$13,a[$4]} }' - ./enh_exp_TStargetgene_exp |sort -u > ./TS_enh_targetgene

#重新处理格式，将TF分割到行，并添加gene表达量和DPM值：
cd /home/ding/all_cmd/script/enh_statistics/
awk 'BEGIN{OFS="\t"}{ split($3,a,","); for(i=1;i<=length(a);i++){print $1,$2,a[i],$4,$5} }' ./TS_enh_TF_targetgene   >./TS_enh_TF_targetgene_1
awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$9"\t"$12]=$10"\t"$13} if(NR>FNR&&($5"\t"$4 in a)){print $0,a[$5"\t"$4]} }' ./enh_exp_TStargetgene_exp ./TS_enh_TF_targetgene_1>./TS_enh_TF_targetgene_2
#再在TS_enh_TF_targetgene_2里加上TF的特异性分数：TF_TSVT文件(TF_TS_score_heatmap.R生成的)提供：将TSPV写入最后一列：
awk 'BEGIN{OFS="\t"}{ if(NR>1){a[$1]=$NF} if(NR>FNR&&$3 in a){print $0,a[$3]} }' ./TF_TSVT ./TS_enh_TF_targetgene_2 >./TS_enh_TF_targetgene_3



###########################################
#获得所有TF与TS enh TS gene的关系数据：
#126个TF与TS Enh的数据来自于：/home/ding/all_cmd/script/enh_statistics/enh_TF/enh_TF_all （enh_exp_TF.sh）

#与enh_exp_targetgene_exp（由./enh_target_gene_exp.sh产生）比较，获得enh调控的特异性的靶基因，(并且靶基因的表达量>0:先不用)
cd /home/ding/all_cmd/script/enh_statistics/

awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1"\t"$3]=$2} if((NR>FNR)&&($10"\t"$13 in a)&&($12!="random")){print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,a[$10"\t"$13]} }' ./pagenbase_ts_gene ./enh_exp_targetgene_exp >./enh_exp_TStargetgene_exp

#生成的enh_exp_TStargetgene_exp 与 all_TF_enh_tissue_1（由enh_exp_TF.sh生成）比较,得到TS enh对应的all TF和对应的TS gene:
cd /home/ding/all_cmd/script/enh_statistics/
awk 'BEGIN{OFS="\t"}{ print $4,$12,$9,$8  }' ./enh_exp_TStargetgene_exp |sort -u | awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1"\t"$2]=$3"\t"$4} if(NR>FNR&&($1"\t"$4 in a)){print $0,a[$1"\t"$4]} }' - ./all_TF_enh_tissue_1 >./TS_enh_TF_targetgene

#特异性enh与特异性的targetgene,(并且target gene的表达量>0,先不用)
#处理上面得到的：./enh_exp_TStargetgene_exp与./enh_spe，结果文件中两列数字分别是gene的表达量和特异性分数。
cd /home/ding/all_cmd/script/enh_statistics/
awk 'BEGIN{OFS="\t"}{ if(NR>1&&$NF=="spe"){print $4,$NF} }' ./enh_spe | awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1]=$2} if(NR>FNR&&$4 in a&&$10>0){print $4,$8,$9,$10,$11,$12,$13,a[$4]} }' - ./enh_exp_TStargetgene_exp |sort -u > ./TS_enh_targetgene


#重新处理格式，将TF分割到行，并添加gene表达量和DPM值：
cd /home/ding/all_cmd/script/enh_statistics/
awk 'BEGIN{OFS="\t"}{ split($3,a,","); for(i=1;i<=length(a);i++){print $1,$2,a[i],$4,$5} }' ./TS_enh_TF_targetgene   >./TS_enh_TF_targetgene_1
awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$9"\t"$12]=$10"\t"$13} if(NR>FNR&&($5"\t"$4 in a)){print $0,a[$5"\t"$4]} }' ./enh_exp_TStargetgene_exp ./TS_enh_TF_targetgene_1>./TS_enh_TF_targetgene_2
#再在TS_enh_TF_targetgene_2里加上TF的特异性分数：TF_TSVT文件(TF_TS_score_heatmap.R生成的)提供：将TSPV写入最后一列：
awk 'BEGIN{OFS="\t"}{ if(NR>1){a[$1]=$NF} if(NR>FNR&&$3 in a){print $0,a[$3]} }' ./TF_TSVT ./TS_enh_TF_targetgene_2 >./TS_enh_TF_targetgene_3














