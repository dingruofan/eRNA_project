#!/bin/bash
#合并各个组织的enh，顺便添加一列作为enh的命名：
#先合并所有的enh包括双向和单向的，命名为enh_x
cd /media/ding/000B49000006264C/eRNA_project/histone
bedops --merge ./Adrenal/enh_find/result.bed ./Brain/enh_find/result.bed ./Breast/enh_find/result.bed ./Heart/enh_find/result.bed ./Liver/enh_find/result.bed ./Lung/enh_find/result.bed ./Ovary/enh_find/result.bed ./Placenta/enh_find/result.bed ./SkeletalMuscle/enh_find/result.bed ./Kidney/enh_find/result.bed | awk 'BEGIN{a=1}{print $0"\tenh_"a;a=a+1}' -  > /home/ding/all_cmd/script/enh_statistics/enh_all_cluster.bed

#获得单项转录本的merge，但因为有的enh在组织1里有bid转录本，但在其他在组织中无bid转录本，这样的enh作为bid_enh，所以需要去掉与bid_enh重合的enh
cd /media/ding/000B49000006264C/eRNA_project/histone
bedops --merge ./Adrenal/enh_find_1/Dnase_histone_unbidir_enh.bed ./Brain/enh_find_1/Dnase_histone_unbidir_enh.bed ./Breast/enh_find_1/Dnase_histone_unbidir_enh.bed ./Heart/enh_find_1/Dnase_histone_unbidir_enh.bed ./Liver/enh_find_1/Dnase_histone_unbidir_enh.bed ./Lung/enh_find_1/Dnase_histone_unbidir_enh.bed ./Ovary/enh_find_1/Dnase_histone_unbidir_enh.bed ./Placenta/enh_find_1/Dnase_histone_unbidir_enh.bed ./SkeletalMuscle/enh_find_1/Dnase_histone_unbidir_enh.bed ./Kidney/enh_find_1/Dnase_histone_unbidir_enh.bed | awk 'BEGIN{}{ if(NR==FNR){a[$1"\t"$2"\t"$3]=$4;next} if(NR>FNR){ b=$1"\t"$2"\t"$3;if(b in a){print b"\t"a[b]}else{print $0"\t无"} }  }' /home/ding/all_cmd/script/enh_statistics/enh_all_cluster.bed -  > /home/ding/all_cmd/script/enh_statistics/enh_unbid_cluster.bed 

#去掉与bid_enh有交集的：
cd /home/ding/all_cmd/script/enh_statistics/
bedops --not-element-of ./enh_unbid_cluster.bed ./enh_bid_cluster.bed > a 
awk '{print $0}' a > ./enh_unbid_cluster.bed


#获得没有eRNA的enh:
cd /home/ding/all_cmd/script/enh_statistics/
bedops --not-element-of ./enh_all_cluster.bed ./enh_bid_cluster.bed | bedops --not-element-of - ./enh_unbid_cluster.bed > ./no_eRNA_enh.bed



#将上面的命名添加到最终的结果中。列说明如下：
#第1-3列为DHS位点也就是enh位点。
#第4列为enh的命名
#第5列为strand
#第6列为双向转录/单向转录位点
#第7-9列为转录因子结合为点
#第10/11列为转录因子名称和表达量

#处理Dnase_histone_unbidir_enh.bed为enh_unbid.bed，也就是给enh添加enh的命名。
#由于既是双向又是单向的enh划分为双向enh，所以单向
cd /media/ding/000B49000006264C/eRNA_project/histone
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  echo $i"中的单向enh聚类"
  awk '{ if(NR==FNR){a[$1"\t"$2"\t"$3]=$4} if(NR>FNR&&($1"\t"$2"\t"$3 in a)){print $1"\t"$2"\t"$3"\t"a[$1"\t"$2"\t"$3]"\t"$5"\t"$4"\t"$6} }' /home/ding/all_cmd/script/enh_statistics/enh_unbid_cluster.bed ./$i/enh_find_1/Dnase_histone_unbidir_enh.bed >./$i/enh_find_1/enh_unbid.bed
done



#处理Dnase_histone_unbidir_enh_TFBS_exp.bed文件为enh_unbid_TF.bed
cd /media/ding/000B49000006264C/eRNA_project/histone
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  awk '{if(NR==FNR){a[$1"\t"$2"\t"$3]=$4} if(NR>FNR&&($1"\t"$2"\t"$3 in a)){print $1"\t"$2"\t"$3"\t"a[$1"\t"$2"\t"$3]"\t"$5"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 } }' /home/ding/all_cmd/script/enh_statistics/enh_unbid_cluster.bed ./$i/enh_find_1/Dnase_histone_unbidir_enh_TFBS_exp.bed >./$i/enh_find_1/enh_unbid_TF.bed

done 



#————————————————————————————————————————————————————————————————————————————————————————————
#统计enh在各个组织中的分布：
cd /media/ding/000B49000006264C/eRNA_project/histone
#首先将所有的enh写入/home/ding/all_cmd/script/enh_statistics/enh_unbid_count
awk 'BEGIN{print "chr\tchr_start\tchr_end\tenhancer"}{print $0}' /home/ding/all_cmd/script/enh_statistics/enh_unbid_cluster.bed > /home/ding/all_cmd/script/enh_statistics/enh_unbid_count

#数组存放要统计的组织：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  echo $i
  #将每一个组织的统计结果写入enh_unbid_count（添加一列）
  #代码解释如下，用了两个awk管道连接，第一个awk是统计这个组织的enh的数目，有记为1,没有记为0.第二个awk将统计好的信息写入文件。
  awk 'BEGIN{}{  if(NR==FNR){a[NR]=$4} if(NR>FNR){b[$4]=$1"\t"$2"\t"$3} }END{ for(i=1;i<=length(a);i++){ if(a[i] in b){print a[i]"\t"1}else{print a[i]"\t"0} } }' /home/ding/all_cmd/script/enh_statistics/enh_unbid_cluster.bed ./$i/enh_find_1/enh_unbid.bed | awk -v tissue="$i"  '{ if(NR==FNR){a[$1]=$2} if(NR>1&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$4]} }' - /home/ding/all_cmd/script/enh_statistics/enh_unbid_count > /home/ding/all_cmd/script/enh_statistics/a
  #由于awk不允许直接覆盖处理的文件，END貌似可以，所以先将结果先保存到临时文件a，然后在复制到统计文件,最后删除a。
  awk '{print $0}' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_unbid_count 
  rm /home/ding/all_cmd/script/enh_statistics/a
done



#统计每个enh在几个组织中存在
awk 'BEGIN{}{ if(NR>1){sum=0;for(i=5;i<=NF;i++){sum=$i+sum} print $1"\t"$2"\t"$3"\t"$4"\t"sum} }END{}' /home/ding/all_cmd/script/enh_statistics/enh_unbid_count |awk '{ if(NR==FNR){a[$4]=$5} if(NR>FNR&&FNR==1){print $0"\ttissue_statistics"} if(NR>FNR&&FNR>1){print $0"\t"a[$4]} }' - /home/ding/all_cmd/script/enh_statistics/enh_unbid_count >/home/ding/all_cmd/script/enh_statistics/a

awk '{print $0}' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_unbid_count 

rm /home/ding/all_cmd/script/enh_statistics/a

gedit /home/ding/all_cmd/script/enh_statistics/enh_unbid_count 

cd /media/ding/000B49000006264C/eRNA_project/histone
awk '{ if(NR==1){print} if($NF==10||$NF==9||$NF==8){print } }' /home/ding/all_cmd/script/enh_statistics/enh_unbid_count 


#————————————————————————————————————————————————————————————————————————————————————————————
#首先统计所有组织，没有eRNA的增强子：
cd /home/ding/all_cmd/script/enh_statistics
bedops --not-element-of  ./enh_all_cluster.bed  ./enh_bid_cluster.bed |bedops --not-element-of - ./enh_unbid_cluster.bed > ./enh_no_eRNA_cluster.bed

#获得各个组织中特异性的no eRNA的增强子：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  bedops --not-element-of ./enh_find/result.bed /home/ding/all_cmd/script/enh_statistics/enh_bid_cluster.bed |bedops --not-element-of - /home/ding/all_cmd/script/enh_statistics/enh_unbid_cluster.bed >./enh_find/enh_no_eRNA.bed
done



#————————————————————————————————————————————————————————————————————————————————————————————
#对各个组织中result.bed和enh_no_eRNA.bed进行名称合并，根据的是enh_all_cluster.bed
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  awk '{ if(NR==FNR){a[$1"\t"$2"\t"$3]=$4;next} if(NR>FNR){print $1"\t"$2"\t"$3"\t"a[$1"\t"$2"\t"$3]} }'  /home/ding/all_cmd/script/enh_statistics/enh_all_cluster.bed   ./enh_find/result.bed > ./enh_find/result_named.bed
   awk '{ if(NR==FNR){a[$1"\t"$2"\t"$3]=$4;next} if(NR>FNR){print $1"\t"$2"\t"$3"\t"a[$1"\t"$2"\t"$3]} }'  /home/ding/all_cmd/script/enh_statistics/enh_all_cluster.bed   ./enh_find/enh_no_eRNA.bed> ./enh_find/enh_no_eRNA_named.bed
done


#————————————————————————————————————————————————————————————————————————————————————————————
#统计no_eRNA_enh在各个组织中的分布：
cd /media/ding/000B49000006264C/eRNA_project/histone
#首先将所有的enh写入/home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_count
awk 'BEGIN{print "chr\tchr_start\tchr_end\tenhancer"}{print $0}' /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_cluster.bed > /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_count


cd /media/ding/000B49000006264C/eRNA_project/histone
#数组存放要统计的组织：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  echo $i
  #将每一个组织的统计结果写入enh_bid_count（添加一列）
  #代码解释如下，用了两个awk管道连接，第一个awk是统计这个组织的enh的数目，有记为1,没有记为0.第二个awk将统计好的信息写入文件。
  awk 'BEGIN{}{  if(NR==FNR){a[NR]=$4} if(NR>FNR){b[$4]=$1"\t"$2"\t"$3} }END{ for(i=1;i<=length(a);i++){ if(a[i] in b){print a[i]"\t"1}else{print a[i]"\t"0} } }' /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_cluster.bed ./$i/enh_find/enh_no_eRNA_named.bed | awk -v tissue="$i"  '{ if(NR==FNR){a[$1]=$2} if(NR>1&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$4]} }' - /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_count >/home/ding/all_cmd/script/enh_statistics/a
  #由于awk不允许直接覆盖处理的文件，END貌似可以，所以先将结果先保存到临时文件a，然后在复制到统计文件,最后删除a。
  awk '{print $0}' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_count 
  rm /home/ding/all_cmd/script/enh_statistics/a
done



#统计每个enh在几个组织中存在
awk 'BEGIN{}{ if(NR>1){sum=0;for(i=5;i<=NF;i++){sum=$i+sum} print $1"\t"$2"\t"$3"\t"$4"\t"sum} }END{}' /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_count |awk '{ if(NR==FNR){a[$4]=$5} if(NR>FNR&&FNR==1){print $0"\ttissue_statistics"} if(NR>FNR&&FNR>1){print $0"\t"a[$4]} }' - /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_count >/home/ding/all_cmd/script/enh_statistics/a

awk '{print $0}' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_count 
rm /home/ding/all_cmd/script/enh_statistics/a

gedit /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_count 

cd /media/ding/000B49000006264C/eRNA_project/histone
awk '{ if(NR==1){print} if($NF==9||$NF==8||$NF==10){print } }' /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_count 




#数组存放要统计的组织：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
#将统计信息写入文件：
cat >>/home/ding/all_cmd/script/enh_statistics/unbid_statistics<<EOF
$i	`awk 'END{print NR}' ./enh_find_1/result.bed`	`awk 'END{print NR}' ./enh_find_1/enh_unbid.bed`	`awk 'END{print NR}' ./enh_find_1/Dnase_histone_unbidir_enh_TFBS_exp.bed` 
EOF


echo -n "符合组蛋白修饰的有："
wc -l ./enh_find_1/result.bed
echo -n "enh共有："
wc -l ./enh_find_1/enh_unbid.bed
echo -n "转录因子为点有："
wc -l ./enh_find_1/Dnase_histone_unbidir_enh_TFBS_exp.bed  

done












