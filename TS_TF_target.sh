#!/bin/bash

######################################
#插入一段寻找每个TS enh 的平均DHS信号的程序：
######################################
#首先获得要查询的tissue和enh的信息，然后从每个组织的DHS中提取其位点信号的均值：
cd /home/ding/all_cmd/script
awk 'BEGIN{OFS="\t"}{ print $1,$4 }' ./enh_statistics/TS_enh_TF_targetgene_3|sort -u |awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$4]=$0} if(NR>FNR){ print a[$1],$2 } }' /home/ding/all_cmd/script/enh_statistics/enh_all_cluster.bed   - >./a

#定义函数获得对应组织的enh的DHS的平均信号：
cd /home/ding/all_cmd/script
#创建空文件用来接收结果：
:>./enh_statistics/TS_enh_TF_targetgene_3_DHS
cat ./a |while read line
do
  #将每个增强子存为一个bed文件，然后使用bwtool进行读取。
  awk -v line1="$line" 'BEGIN{OFS="\t";split(line1,a,"\t"); print a[1],a[2],a[3],a[4]}' >./a_enh.bed
  tissue=`awk -v tissue="$line" 'BEGIN{OFS="\t";split(tissue,a,"\t");print a[5]}' `
  #进入要操作的目录：
  DHS_mean_enh=`bwtool summary ./a_enh.bed /media/ding/000B49000006264C/eRNA_project/DNase/$tissue/new_DHS_TPM.bw  /dev/stdout -with-sum -fill=0|awk '{print $8}' -`
  awk -v line1="$line" -v DHS_mean_enh="$DHS_mean_enh" 'BEGIN{OFS="\t";print line1,DHS_mean_enh}' >>./enh_statistics/TS_enh_TF_targetgene_3_DHS
  echo "已经完成"$line"	$DHS_mean_enh"
done

rm ./a ./a_enh.bed



#################################

cd /home/ding/all_cmd/script

#由于有些TF的名称和基因的名称相同，第三列和第五列。在相同名称的基因后加"Gene-"
#说明：先找到相同的gene-TF名再替换：
awk 'BEGIN{OFS='\t'}{ a[$3];b[$5] }END{ for(i in a){if(i in b){print i}} }' ./enh_statistics/TS_enh_TF_targetgene_3 | awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1]} if(NR>FNR){ if($5 in a){print $1,$2,$3,$4,$5"_",$6,$7,$8}else{print $0} } }' - ./enh_statistics/TS_enh_TF_targetgene_3 > ./enh_statistics/TS_enh_TF_targetgene_3.1


#将DHS信号加入./enh_statistics/TS_enh_TF_targetgene_3.1的后面并按照最后一列DHS排序，生成3.2
awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$4"\t"$5]=$NF} if(NR>FNR){print $0,a[$1"\t"$4]} }' ./enh_statistics/TS_enh_TF_targetgene_3_DHS ./enh_statistics/TS_enh_TF_targetgene_3.1 |sort -nk 9 > ./enh_statistics/TS_enh_TF_targetgene_3.2




############
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
for i in ${tissue[@]}
do
  #先创建空文件：
  :>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_relation"

  #添加 TF	target	TSVP	"gene"
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if($4==tissue){print $3,$5,"TF","gene","TF-gene"} }' ./enh_statistics/TS_enh_TF_targetgene_3.2 |sort -u - >./enh_statistics/cytoscape/a
  awk -v tissue=$i 'BEGIN{OFS="\t";print "start","end","start_label","end_label","relation"}{ if(NR==FNR){a[$1"\t"$2]} if(NR>FNR){ split($2,b,"_"); if($1"\t"b[1] in a||b[1]"\t"$1 in a)print $0 }}' /media/ding/000B49000006264C/eRNA_project/biogrid/pro_relation ./enh_statistics/cytoscape/a >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_relation"
  rm ./enh_statistics/cytoscape/a

  #添加 enh	TF	"enhancer"	TSVP
  #print "enhancer","TF","enh_label","TF_label","relation"
  awk -v tissue=$i 'BEGIN{OFS="\t";}{ if($4==tissue){print $1,$3,"enhancer","TF","enhancer-TF"} }' ./enh_statistics/TS_enh_TF_targetgene_3.2 |sort -u - >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_relation"
  
  #添加 enh	target	"enhancer"	"gene"
  #print "enhancer","gene","enh_label","gene_label","relation"
  awk -v tissue=$i 'BEGIN{OFS="\t";}{ if($4==tissue){print $1,$5,"enhancer","gene","enhancer-gene"} }' ./enh_statistics/TS_enh_TF_targetgene_3.2 |sort -u - >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_relation"

  echo $i
done




#对value注释，用于划分形状和颜色：
cd /home/ding/all_cmd/script
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
for i in ${tissue[@]}
do
  #先创建空文件：
  :>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_value"
  awk 'BEGIN{OFS="\t";print "new_start","value"}'  >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_value"
  #添加TF	TSVP
  awk -v tissue=$i 'BEGIN{OFS="\t";}{ if($4==tissue){print $3,$8} }' ./enh_statistics/TS_enh_TF_targetgene_3.2 |sort -u - >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_value"
  #添加gene	100
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if($4==tissue){print $5,100} }' ./enh_statistics/TS_enh_TF_targetgene_3.2 |sort -u - >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_value"

  #enhancer	信号平均值：来源于上面生成的./enh_statistics/TS_enh_TF_targetgene_3_DHS，表示增强子活性：
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if($4==tissue){print $1,$9} }' ./enh_statistics/TS_enh_TF_targetgene_3.2 |sort -u - >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_value"
  #awk -v tissue="$i" 'BEGIN{OFS="\t"}{ if($5==tissue)print $4,$6 }' ./enh_statistics/TS_enh_TF_targetgene_3_DHS >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_value"
  echo $i
done





#替换每个组织的TSPV为TSVT：./enh_statistics/TF_TSVT
cd /home/ding/all_cmd/script
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
for i in ${tissue[@]}
do
  awk -v tissue=$i 'BEGIN{OFS="\t";tissue_id=12}{if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}}  if(NR>1&&NR==FNR){a[$1]=$tissue_id} if(NR>FNR&&FNR==1){print $0} if(NR>FNR&&FNR>1&&$1 in a&&$2<0){print $1,a[$1]} if(NR>FNR&&FNR>1&&$!(1 in a)&&!($2<0)){print $0} }' ./enh_statistics/TF_TSVT /media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_value" >./enh_statistics/a
  awk '{print $0}' ./enh_statistics/a > /media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_value"
  rm ./enh_statistics/a
  echo $i
done





#目的：cytoscape圆圈排序，实现增强子和gene一一对应的关系。
#获得每个节点的位置，enh和gene的位置相同。
cd /home/ding/all_cmd/script
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
for i in ${tissue[@]}
do
  #先创建空文件：
  :>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_location"
  awk 'BEGIN{OFS="\t";print "node","location"}'  >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_location"
  #添加TF	NR
  awk -v tissue=$i 'BEGIN{OFS="\t";}{ if($4==tissue){print $3,NR} }' ./enh_statistics/TS_enh_TF_targetgene_3.2 |sort -u - >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_location"
  #添加gene	NR
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if($4==tissue){print $5,NR} }' ./enh_statistics/TS_enh_TF_targetgene_3.2 |sort -u - >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_location"
  #enhancer	NR
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if($4==tissue){print $1,NR} }' ./enh_statistics/TS_enh_TF_targetgene_3.2 |sort -u - >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_location"
  echo $i
done




#目的：节点上的edge越多，其label size或者自体越大。
#获得每个节点的edge的数目。
cd /home/ding/all_cmd/script
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
#tissue=(Adrenal)
for i in ${tissue[@]}
do
  #先创建空文件：
  :>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_edge_num"
  awk 'BEGIN{OFS="\t";print "node","edge_num"}'  >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_edge_num"
  #TF edge_num
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if($4==tissue){ a[$3]++; } }END{ for(i in a){print i,a[i]} }' ./enh_statistics/TS_enh_TF_targetgene_3.2  >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_edge_num"
  #gene edge_num
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if($4==tissue){ a[$5]++; } }END{ for(i in a){print i,a[i]} }' ./enh_statistics/TS_enh_TF_targetgene_3.2 >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_edge_num"
  #enh edge_num
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if($4==tissue){ a[$1]++; } }END{ for(i in a){print i,a[i]} }' ./enh_statistics/TS_enh_TF_targetgene_3.2 >>/media/ding/000B49000006264C/eRNA_project/cytoscape/$i"_enh_TF_target_edge_num"
  echo $i
done


















