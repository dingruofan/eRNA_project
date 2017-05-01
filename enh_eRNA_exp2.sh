#先将ctss.bed转换为标准的bed6格式，然后使用bedtools中的bedToBam转为bam文件,然后再使用Rseqc的RPKM_count来计算每个组织：
#对TCs进行排序并生成bed12格式的TCs，只有bed12才可以使用Rseqc程序
cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
sort-bed ./cage_enh/TCs.bed>./cage_enh/TCs_sorted.bed
awk '{ if($1=="chrM"){next}else{ print $0"\t"$2"\t"$3"\t0,0,0\t1\t"$3-$2"\t0"} }' ./cage_enh/TCs_sorted.bed>./cage_enh/TCs_sorted.bed12


tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  cd /media/ding/000B49000006264C/eRNA_project/fantom/$i"_CAGE"
  #获得所有的ctss：
  awk 'BEGIN{OFS="\t"} { a[NR]=$5;b[NR]=$0; }END{ for(i=1;i<=length(a);i++){ for(j=1;j<=a[i];j++){print b[i]} } }' ./ctss.bed >./ctss_all.bed
  #转为bam文件：
  bedToBam -i ./ctss_all.bed -g ../../hg19.chrom.sizes >./ctss_all.bam
  samtools sort ./ctss_all.bam ./ctss_all_sorted
  samtools index ./ctss_all_sorted.bam ./ctss_all_sorted.bam.bai
  #获得所有双向转录本分割为单个转录本的bed12格式文件：
  cd /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE
  bed12ToBed6 -i ./cage_enh/bidir.pairs.bed | awk 'BEGIN{OFS="\t";} { if(NR%2==1){print $1,$2,$3,$4",-",$5,"-"} if(NR%2==0){print $1,$2,$3,$4",+",$5,"+"} } ' - |sortBed -i - >./cage_enh/bidir_pairs_sorted.bed6
  #转为bed12格式：
  awk '{ if($1=="chrM"){next}else{ print $0"\t"$2"\t"$3"\t0,0,0\t1\t"$3-$2"\t0"} }' ./cage_enh/bidir_pairs_sorted.bed6 >./cage_enh/bidir_pairs_sorted.bed12
  #使用RSeqc的工具获得每个组织在每个转录本的RPKM：
  cd /media/ding/000B49000006264C/eRNA_project/fantom/$i"_CAGE"
  /home/ding/tool/RSeQC-2.3.7/scripts/RPKM_count.py -d '++,–' -r /media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/cage_enh/TCs.bed12  -i  ./ctss_all_sorted.bam -o bid_trans_RPKM

done


###########################
#已经获得了/media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/cage_enh/bidir_pairs_expression_TPM.matrix双向每个转录本的表达量
#单向的：/media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/cage_enh/TCs_sorted_exp_tpm.matrix的转录本的表达量

  #创建/home/ding/all_cmd/script/enh_statistics/enh_bid_cluster_for_exp.bed
  awk 'BEGIN{OFS="\t";}{ print $1,$2,$3,$4",-";print $1,$2,$3,$4",+";}' /home/ding/all_cmd/script/enh_statistics/enh_bid_cluster.bed >/home/ding/all_cmd/script/enh_statistics/enh_bid_cluster_for_exp.bed
  awk 'BEGIN{print "chr\tchr_start\tchr_end\tenhancer"}{print $0}' /home/ding/all_cmd/script/enh_statistics/enh_bid_cluster_for_exp.bed > /home/ding/all_cmd/script/enh_statistics/enh_bid_exp2

tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  #首先获得每个组织的eRNA的表达量放在各个组织文件夹下，注意同一个enh的在不同组织的eRNA可能不同。
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  awk -v tissue=$i 'BEGIN{tissue_ID=0;OFS="\t"}{ if(NR==1){ for(i=1;i<=NF;i++){ if($i==tissue){tissue_ID=i} } } if((NR==FNR)&&(NR>1)){a[$1]=$tissue_ID}  if(NR>FNR&&($6",-" in a)){print $1,$2,$3,$4,$5,"-",$7,a[$6",-"]} if(NR>FNR&&($6",+" in a)){print $1,$2,$3,$4,$5,"+",$7,a[$6",+"]}}'  ../../fantom/All_tissue_CAGE/bid_cage_enh_quantify/bidir_pairs_expression_tpm.matrix  ./enh_find/enh_bid.bed >./enh_find/enh_bid_exp2.bed 


 #然后再将表达量合并到enh_bid_exp2，没有的表达量为0。
  awk -v tissue=$i '{   if(NR==FNR){ a[$4","$6]=$NF }  if(NR>FNR){ if($4 in a){ print  $1,$2,$3,$4,a[$4] }else{print  $1,$2,$3,$4,0} }   }' ./enh_find/enh_bid_exp2.bed   /home/ding/all_cmd/script/enh_statistics/enh_bid_cluster_for_exp.bed |awk  -v tissue=$i 'BEGIN{OFS="\t"}{  if(NR==FNR){a[$1"\t"$2"\t"$3"\t"$4]=$5} if(NR>FNR&&FNR==1){print $0"\t"tissue}  if(NR>FNR&&FNR>1){print $0,a[$1"\t"$2"\t"$3"\t"$4]} }' -  /home/ding/all_cmd/script/enh_statistics/enh_bid_exp2 >/home/ding/all_cmd/script/enh_statistics/a 
  awk '{ print $0 }' /home/ding/all_cmd/script/enh_statistics/a > /home/ding/all_cmd/script/enh_statistics/enh_bid_exp2
  rm /home/ding/all_cmd/script/enh_statistics/a
  echo $i
done























#————————————————————————————————————————————————————————————————————————————————————————————————————————————————————


#做unbid eRAN的表达量分析：
#从所有的enh_unbid_cluster.bed匹配到每个组织的eRNA，然后匹配到表达矩阵。








