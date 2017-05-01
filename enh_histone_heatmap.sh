#!/bin/bash

#######################################################################
#各个组织bid/unbid/no_eRNA增强子中心附近2kb的组蛋白修饰热图：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  #获得每个组织的bid/unbid/no_eRNA增强子中心附近2kb附近3种组蛋白修饰的平均信号,貌似bwtool不支持管道。。
  

 #读入该组织总的信号值：
function read_H3K4me1_sigsum(){
 cat ./bg/H3K4me1_sigsum |while read line; do 
   echo `awk 'BEGIN{ print '"$line"' }'` 
   done
}
function read_H3K4me3_sigsum(){
 cat ./bg/H3K4me3_sigsum |while read line; do 
   echo `awk 'BEGIN{ print '"$line"' }'` 
   done
}
function read_H3K27ac_sigsum(){
 cat ./bg/H3K27ac_sigsum |while read line; do 
   echo `awk 'BEGIN{ print '"$line"' }'` 
   done
}
H3K4me1_sig=`read_H3K4me1_sigsum`
H3K4me3_sig=`read_H3K4me3_sigsum`
H3K27ac_sig=`read_H3K27ac_sigsum`

echo $H3K4me1_sig
echo $H3K4me3_sig
echo $H3K27ac_sig
 #处理每个组织的bid信号:
  awk '{ print $1"\t"$2"\t"$3 }'  ./enh_find/enh_bid.bed >./enh_find/a
  bwtool agg 2000:2000  ./enh_find/a  ./bg/H3K4me1.bw,./bg/H3K4me3.bw,./bg/H3K27ac.bw  /dev/stdout > ./enh_find/enh_bid_histone_signal
 #进行归一化：除以总的信号数
 awk -v H3K4me1_sig="$H3K4me1_sig" -v H3K4me3_sig="$H3K4me3_sig" -v H3K27ac_sig="$H3K27ac_sig" '{ printf("%d\t%0.9f\t%0.9f\t%0.9f\n",$1,$2*10^6/H3K4me1_sig,$3*10^6/H3K4me3_sig,$4*10^6/H3K27ac_sig) }' ./enh_find/enh_bid_histone_signal >./enh_find/a
 awk '{print $0}' ./enh_find/a > ./enh_find/enh_bid_histone_signal


 #处理每个组织的no_eRNA信号:
  awk '{ print $1"\t"$2"\t"$3 }'  ./enh_find/enh_no_eRNA_named.bed>./enh_find/a
  bwtool agg 2000:2000  ./enh_find/a  ./bg/H3K4me1.bw,./bg/H3K4me3.bw,./bg/H3K27ac.bw  /dev/stdout > ./enh_find/enh_no_eRNA_histone_signal
 #进行归一化：除以总的信号数
 awk -v H3K4me1_sig="$H3K4me1_sig" -v H3K4me3_sig="$H3K4me3_sig" -v H3K27ac_sig="$H3K27ac_sig" '{ printf("%d\t%0.9f\t%0.9f\t%0.9f\n",$1,$2*10^6/H3K4me1_sig,$3*10^6/H3K4me3_sig,$4*10^6/H3K27ac_sig) }'  ./enh_find/enh_no_eRNA_histone_signal >./enh_find/a
 awk '{print $0}' ./enh_find/a > ./enh_find/enh_no_eRNA_histone_signal


 #处理每个组织的unbid信号:
  awk '{ print $1"\t"$2"\t"$3 }'  ./enh_find_1/enh_unbid.bed>./enh_find/a
  bwtool agg 2000:2000  ./enh_find/a  ./bg/H3K4me1.bw,./bg/H3K4me3.bw,./bg/H3K27ac.bw  /dev/stdout > ./enh_find_1/enh_unbid_histone_signal
 #进行归一化：除以总的信号数
 awk -v H3K4me1_sig="$H3K4me1_sig" -v H3K4me3_sig="$H3K4me3_sig" -v H3K27ac_sig="$H3K27ac_sig" '{ printf("%d\t%0.9f\t%0.9f\t%0.9f\n",$1,$2*10^6/H3K4me1_sig,$3*10^6/H3K4me3_sig,$4*10^6/H3K27ac_sig) }' ./enh_find_1/enh_unbid_histone_signal>./enh_find/a
 awk '{print $0}' ./enh_find/a > ./enh_find_1/enh_unbid_histone_signal
 echo $i
done



#合并各个组织的组蛋白修饰信号到一个文件，按照组蛋白修饰：
#######################################################################


#生成bid信号的三个文件
 awk 'BEGIN{print "location"} { print $1 }' /media/ding/000B49000006264C/eRNA_project/histone/Adrenal/enh_find/enh_bid_histone_signal |tee /home/ding/all_cmd/script/enh_statistics/enh_bid_H3K4me1_signal   /home/ding/all_cmd/script/enh_statistics/enh_bid_H3K4me3_signal  /home/ding/all_cmd/script/enh_statistics/enh_bid_H3K27ac_signal 
#生成unbid信号的三个文件
 awk 'BEGIN{print "location"} { print $1 }' /media/ding/000B49000006264C/eRNA_project/histone/Adrenal/enh_find_1/enh_unbid_histone_signal |tee /home/ding/all_cmd/script/enh_statistics/enh_unbid_H3K4me1_signal   /home/ding/all_cmd/script/enh_statistics/enh_unbid_H3K4me3_signal  /home/ding/all_cmd/script/enh_statistics/enh_unbid_H3K27ac_signal 
#生成no_eRNA信号的三个文件
 awk 'BEGIN{print "location"} { print $1 }' /media/ding/000B49000006264C/eRNA_project/histone/Adrenal/enh_find/enh_no_eRNA_histone_signal |tee /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_H3K4me1_signal   /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_H3K4me3_signal  /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_H3K27ac_signal 

##################
#
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  # bid H3K4me1:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$2 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find/enh_bid_histone_signal  /home/ding/all_cmd/script/enh_statistics/enh_bid_H3K4me1_signal >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_bid_H3K4me1_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  # bid  H3K4me3:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$3 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find/enh_bid_histone_signal  /home/ding/all_cmd/script/enh_statistics/enh_bid_H3K4me3_signal >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_bid_H3K4me3_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  #  bid H3K27ac:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$4 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find/enh_bid_histone_signal  /home/ding/all_cmd/script/enh_statistics//enh_bid_H3K27ac_signal  >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_bid_H3K27ac_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
done

tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  # unbid H3K4me1:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$2 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find_1/enh_unbid_histone_signal  /home/ding/all_cmd/script/enh_statistics/enh_unbid_H3K4me1_signal >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_unbid_H3K4me1_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  # unbid  H3K4me3:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$3 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find_1/enh_unbid_histone_signal  /home/ding/all_cmd/script/enh_statistics/enh_unbid_H3K4me3_signal >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_unbid_H3K4me3_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  #  unbid H3K27ac:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$4 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find_1/enh_unbid_histone_signal  /home/ding/all_cmd/script/enh_statistics//enh_unbid_H3K27ac_signal  >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_unbid_H3K27ac_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
done


tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  # no_eRNA H3K4me1:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$2 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find/enh_no_eRNA_histone_signal  /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_H3K4me1_signal >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_H3K4me1_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  # no_eRNA  H3K4me3:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$3 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find/enh_no_eRNA_histone_signal  /home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_H3K4me3_signal >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_H3K4me3_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  #  no_eRNA H3K27ac:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$4 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find/enh_no_eRNA_histone_signal  /home/ding/all_cmd/script/enh_statistics//enh_no_eRNA_H3K27ac_signal  >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_no_eRNA_H3K27ac_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a

echo $i
done


 #合并bid各个组织的的组蛋白修饰信号到一个文件：
cd /home/ding/all_cmd/script/enh_statistics
awk 'BEGIN{OFS="\t"} (ARGIND==1){a[$1]=$0;}  (ARGIND==2){b[$1]=$0} (ARGIND==3){c[$1]=$0} END{  range_start=(length(a)-1)/2*-1; range_end=(length(a)-1)/2;   {print a["location"],"type" }  for(i=range_start;i<=range_end;i++){if(i==0)continue;print a[i],"bid"}  for(i=range_start;i<=range_end;i++){if(i==0)continue;print b[i],"unbid"}  for(i=range_start;i<=range_end;i++){if(i==0)continue;print c[i],"no_eRNA"} } ' ./enh_bid_H3K4me1_signal  ./enh_unbid_H3K4me1_signal  ./enh_no_eRNA_H3K4me1_signal >./enh_H3K4me1_signal

awk 'BEGIN{OFS="\t"} (ARGIND==1){a[$1]=$0;}  (ARGIND==2){b[$1]=$0} (ARGIND==3){c[$1]=$0} END{  range_start=(length(a)-1)/2*-1; range_end=(length(a)-1)/2;   {print a["location"],"type" }  for(i=range_start;i<=range_end;i++){if(i==0)continue;print a[i],"bid"}  for(i=range_start;i<=range_end;i++){if(i==0)continue;print b[i],"unbid"}  for(i=range_start;i<=range_end;i++){if(i==0)continue;print c[i],"no_eRNA"} } ' ./enh_bid_H3K4me3_signal  ./enh_unbid_H3K4me3_signal  ./enh_no_eRNA_H3K4me3_signal >./enh_H3K4me3_signal

awk 'BEGIN{OFS="\t"} (ARGIND==1){a[$1]=$0;}  (ARGIND==2){b[$1]=$0} (ARGIND==3){c[$1]=$0} END{  range_start=(length(a)-1)/2*-1; range_end=(length(a)-1)/2;   {print a["location"],"type" }  for(i=range_start;i<=range_end;i++){if(i==0)continue;print a[i],"bid"}  for(i=range_start;i<=range_end;i++){if(i==0)continue;print b[i],"unbid"}  for(i=range_start;i<=range_end;i++){if(i==0)continue;print c[i],"no_eRNA"} } ' ./enh_bid_H3K27ac_signal  ./enh_unbid_H3K27ac_signal  ./enh_no_eRNA_H3K27ac_signal >./enh_H3K27ac_signal

rm ./enh_bid_H3K4me1_signal  ./enh_unbid_H3K4me1_signal  ./enh_no_eRNA_H3K4me1_signal ./enh_bid_H3K4me3_signal  ./enh_unbid_H3K4me3_signal  ./enh_no_eRNA_H3K4me3_signal  ./enh_bid_H3K27ac_signal  ./enh_unbid_H3K27ac_signal  ./enh_no_eRNA_H3K27ac_signal 


###########################################
#生成wig格式的文件倒入igv即可。
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  #各个组织H3K4me1 bid的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=255,0,0\nfixedStep chrom=chr1 start=1 setp=1 span=1"} { if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if(NR>1&&$12=="bid")print $tissue_id }'   enh_H3K4me1_signal >./histone_wig/H3K4me1/bid/$i.wig
  #各个组织H3K4me1 unbid的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=255,0,0\nfixedStep chrom=chr1 start=1 setp=1 span=1"} {  if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="unbid")print $tissue_id }'   enh_H3K4me1_signal >./histone_wig/H3K4me1/unbid/$i.wig
  #各个组织H3K4me1 no_eRNA的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=255,0,0\nfixedStep chrom=chr1 start=1 setp=1 span=1"} {  if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="no_eRNA")print $tissue_id }'   enh_H3K4me1_signal >./histone_wig/H3K4me1/no_eRNA/$i.wig
##########
  #各个组织H3K4me3 bid的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=0,0,255\nfixedStep chrom=chr1 start=1 setp=1 span=1"} {  if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="bid")print $tissue_id }'   enh_H3K4me3_signal >./histone_wig/H3K4me3/bid/$i.wig
  #各个组织H3K4me3 unbid的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=0,0,255\nfixedStep chrom=chr1 start=1 setp=1 span=1"} {  if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="unbid")print $tissue_id }'   enh_H3K4me3_signal >./histone_wig/H3K4me3/unbid/$i.wig
  #各个组织H3K4me3 no_eRNA的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=0,0,255\nfixedStep chrom=chr1 start=1 setp=1 span=1"} {  if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="no_eRNA")print $tissue_id }'   enh_H3K4me3_signal >./histone_wig/H3K4me3/no_eRNA/$i.wig
##########
  #各个组织H3K27ac bid的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=0,255,0\nfixedStep chrom=chr1 start=1 setp=1 span=1"} { if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="bid")print $tissue_id }'   enh_H3K27ac_signal >./histone_wig/H3K27ac/bid/$i.wig
  #各个组织H3K27ac unbid的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=0,255,0\nfixedStep chrom=chr1 start=1 setp=1 span=1"} {  if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="unbid")print $tissue_id }'   enh_H3K27ac_signal >./histone_wig/H3K27ac/unbid/$i.wig
  #各个组织H3K27ac no_eRNA的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=0,255,0\nfixedStep chrom=chr1 start=1 setp=1 span=1"} {  if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="no_eRNA")print $tissue_id }'   enh_H3K27ac_signal >./histone_wig/H3K27ac/no_eRNA/$i.wig

echo $i
done













