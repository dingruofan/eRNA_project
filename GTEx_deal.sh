cd /media/ding/000B49000006264C/eRNA_project/GTEx

#输出文件到./GTEx.rpkm
#注意：获得单向转录的enh要去除已知的双向转录的enh，一个DHS既有bid CAGE,也有unbid的TCs，算为bid_enh.

awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==3){print "gene_id","gene_name"} if(NR>3){print $1,$2} }' ./GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct >./GTEx_rpkm

#数组存放9个组织,注意只有9个组织：
tissue=(GTEX-1399S-0426-SM-5IFG5 GTEX-1117F-3226-SM-5N9CT GTEX-11P81-1926-SM-5BC53 GTEX-11ZUS-0226-SM-5FQT8 GTEX-11OF3-1326-SM-5N9FJ GTEX-14AS3-0126-SM-5Q5F4 GTEX-11P81-0226-SM-5HL5M GTEX-11P81-1526-SM-5P9GS GTEX-XOT4-0526-SM-4B66O)
#tissue=(Adrenal)
for i in ${tissue[@]}
do 
  awk -v tissue_id=$i 'BEGIN{FS="\t";OFS="\t";j=1}{ if(NR==FNR&&NR==3){for(i=1;i<=NF;i++){if($i==tissue_id){j=i}}}  if(NR==FNR&&NR>3){a[$1]=$j}  if(NR>FNR&&FNR==1){print $0,tissue_id} if(NR>FNR&&FNR>1){print $0,a[$1]} }' ./GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct ./GTEx_rpkm >./a
  cp ./a ./GTEx_rpkm
  rm ./a
  echo $i"完成1"
done


#将上面生成文件的tissue_id替换为tissue:
awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR){a[$2]=$1} if(NR>FNR&&FNR==1){print $1,$2,a[$3],a[$4],a[$5],a[$6],a[$7],a[$8],a[$9],a[$10],a[$11]} if(NR>FNR&&FNR>1){print $0} }' ./tissue_id ./GTEx_rpkm >./a
  cp ./a ./GTEx_rpkm
  rm ./a
  echo "ID转换完成。"

















