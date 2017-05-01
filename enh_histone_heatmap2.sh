#!/bin/bash

#######################################################################
#各个组织spe/uni/other增强子中心附近2kb的组蛋白修饰热图：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  #获得每个组织的spe/uni/other增强子中心附近2kb附近3种组蛋白修饰的平均信号,貌似bwtool不支持管道。。
  
 #处理每个组织的spe信号:
  awk '{ print $1"\t"$2"\t"$3"\t"$4 }'  ./enh_find/enh_bid.bed ./enh_find/enh_no_eRNA_named.bed  ./enh_find_1/enh_unbid.bed | awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR&&NR>1){a[$4]=$7} if(NR>FNR){print $1,$2,$3,$4,a[$4]} }' /home/ding/all_cmd/script/enh_statistics/enh_spe - |awk 'BEGIN{FS="\t";OFS="\t"}{ if($5=="spe"){print $1,$2,$3} }' - >./enh_find/a
  bwtool agg 2000:2000  ./enh_find/a  ./bg/H3K4me1.bw,./bg/H3K4me3.bw,./bg/H3K27ac.bw  /dev/stdout > ./enh_find/enh_spe_histone_signal


 #处理每个组织的uni信号:
    awk '{ print $1"\t"$2"\t"$3"\t"$4 }'  ./enh_find/enh_bid.bed ./enh_find/enh_no_eRNA_named.bed  ./enh_find_1/enh_unbid.bed | awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR&&NR>1){a[$4]=$7} if(NR>FNR){print $1,$2,$3,$4,a[$4]} }' /home/ding/all_cmd/script/enh_statistics/enh_spe - |awk 'BEGIN{FS="\t";OFS="\t"}{ if($5=="uni"){print $1,$2,$3} }' - >./enh_find/a
  bwtool agg 2000:2000  ./enh_find/a  ./bg/H3K4me1.bw,./bg/H3K4me3.bw,./bg/H3K27ac.bw  /dev/stdout > ./enh_find/enh_uni_histone_signal


 #处理每个组织的other信号:
    awk '{ print $1"\t"$2"\t"$3"\t"$4 }'  ./enh_find/enh_bid.bed ./enh_find/enh_no_eRNA_named.bed  ./enh_find_1/enh_unbid.bed | awk 'BEGIN{FS="\t";OFS="\t"}{ if(NR==FNR&&NR>1){a[$4]=$7} if(NR>FNR){print $1,$2,$3,$4,a[$4]} }' /home/ding/all_cmd/script/enh_statistics/enh_spe - |awk 'BEGIN{FS="\t";OFS="\t"}{ if($5=="other"){print $1,$2,$3} }' - >./enh_find/a
  bwtool agg 2000:2000  ./enh_find/a  ./bg/H3K4me1.bw,./bg/H3K4me3.bw,./bg/H3K27ac.bw  /dev/stdout > ./enh_find_1/enh_other_histone_signal
 echo $i
done



#合并各个组织的组蛋白修饰信号到一个文件，按照组蛋白修饰：
#######################################################################


#生成spe信号的三个文件
 awk 'BEGIN{print "location"} { print $1 }' /media/ding/000B49000006264C/eRNA_project/histone/Adrenal/enh_find/enh_spe_histone_signal |tee /home/ding/all_cmd/script/enh_statistics/enh_spe_H3K4me1_signal   /home/ding/all_cmd/script/enh_statistics/enh_spe_H3K4me3_signal  /home/ding/all_cmd/script/enh_statistics/enh_spe_H3K27ac_signal 
#生成uni信号的三个文件
 awk 'BEGIN{print "location"} { print $1 }' /media/ding/000B49000006264C/eRNA_project/histone/Adrenal/enh_find/enh_uni_histone_signal |tee /home/ding/all_cmd/script/enh_statistics/enh_uni_H3K4me1_signal   /home/ding/all_cmd/script/enh_statistics/enh_uni_H3K4me3_signal  /home/ding/all_cmd/script/enh_statistics/enh_uni_H3K27ac_signal 
#生成other信号的三个文件
 awk 'BEGIN{print "location"} { print $1 }' /media/ding/000B49000006264C/eRNA_project/histone/Adrenal/enh_find_1/enh_other_histone_signal |tee /home/ding/all_cmd/script/enh_statistics/enh_other_H3K4me1_signal   /home/ding/all_cmd/script/enh_statistics/enh_other_H3K4me3_signal  /home/ding/all_cmd/script/enh_statistics/enh_other_H3K27ac_signal 

##################
#spe的各个组织的组蛋白修饰：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  # spe H3K4me1:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$2 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find/enh_spe_histone_signal  /home/ding/all_cmd/script/enh_statistics/enh_spe_H3K4me1_signal >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_spe_H3K4me1_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  # spe  H3K4me3:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$3 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find/enh_spe_histone_signal  /home/ding/all_cmd/script/enh_statistics/enh_spe_H3K4me3_signal >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_spe_H3K4me3_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  #  spe H3K27ac:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$4 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find/enh_spe_histone_signal  /home/ding/all_cmd/script/enh_statistics//enh_spe_H3K27ac_signal  >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_spe_H3K27ac_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  echo $i
done

#uni的各个组织的组蛋白修饰：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  # uni H3K4me1:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$2 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find/enh_uni_histone_signal  /home/ding/all_cmd/script/enh_statistics/enh_uni_H3K4me1_signal >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_uni_H3K4me1_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  # uni  H3K4me3:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$3 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find/enh_uni_histone_signal  /home/ding/all_cmd/script/enh_statistics/enh_uni_H3K4me3_signal >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_uni_H3K4me3_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  #  uni H3K27ac:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$4 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find/enh_uni_histone_signal  /home/ding/all_cmd/script/enh_statistics//enh_uni_H3K27ac_signal  >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_uni_H3K27ac_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  echo $i
done

#other的各个组织的组蛋白修饰：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
cd /media/ding/000B49000006264C/eRNA_project/histone/$i
  # other H3K4me1:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$2 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find_1/enh_other_histone_signal  /home/ding/all_cmd/script/enh_statistics/enh_other_H3K4me1_signal >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_other_H3K4me1_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  # other  H3K4me3:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$3 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find_1/enh_other_histone_signal  /home/ding/all_cmd/script/enh_statistics/enh_other_H3K4me3_signal >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_other_H3K4me3_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  #  other H3K27ac:
  awk -v tissue="$i" '{ if(NR==FNR){ a[$1]=$4 } if(NR>FNR&&FNR==1){print $0"\t"tissue} if(NR>FNR&&FNR>1){print $0"\t"a[$1]} }'  ./enh_find_1/enh_other_histone_signal  /home/ding/all_cmd/script/enh_statistics//enh_other_H3K27ac_signal  >/home/ding/all_cmd/script/enh_statistics/a
  awk ' { print $0 } ' /home/ding/all_cmd/script/enh_statistics/a >/home/ding/all_cmd/script/enh_statistics/enh_other_H3K27ac_signal 
  rm /home/ding/all_cmd/script/enh_statistics/a
  echo $i
done


 #合并spe各个组织的的组蛋白修饰信号到一个文件：
cd /home/ding/all_cmd/script/enh_statistics
awk 'BEGIN{OFS="\t"} (ARGIND==1){a[$1]=$0;}  (ARGIND==2){b[$1]=$0} (ARGIND==3){c[$1]=$0} END{  range_start=(length(a)-1)/2*-1; range_end=(length(a)-1)/2;   {print a["location"],"type" }  for(i=range_start;i<=range_end;i++){if(i==0)continue;print a[i],"spe"}  for(i=range_start;i<=range_end;i++){if(i==0)continue;print b[i],"uni"}  for(i=range_start;i<=range_end;i++){if(i==0)continue;print c[i],"other"} } ' ./enh_spe_H3K4me1_signal  ./enh_uni_H3K4me1_signal  ./enh_other_H3K4me1_signal >./enh_H3K4me1_signal

awk 'BEGIN{OFS="\t"} (ARGIND==1){a[$1]=$0;}  (ARGIND==2){b[$1]=$0} (ARGIND==3){c[$1]=$0} END{  range_start=(length(a)-1)/2*-1; range_end=(length(a)-1)/2;   {print a["location"],"type" }  for(i=range_start;i<=range_end;i++){if(i==0)continue;print a[i],"spe"}  for(i=range_start;i<=range_end;i++){if(i==0)continue;print b[i],"uni"}  for(i=range_start;i<=range_end;i++){if(i==0)continue;print c[i],"other"} } ' ./enh_spe_H3K4me3_signal  ./enh_uni_H3K4me3_signal  ./enh_other_H3K4me3_signal >./enh_H3K4me3_signal

awk 'BEGIN{OFS="\t"} (ARGIND==1){a[$1]=$0;}  (ARGIND==2){b[$1]=$0} (ARGIND==3){c[$1]=$0} END{  range_start=(length(a)-1)/2*-1; range_end=(length(a)-1)/2;   {print a["location"],"type" }  for(i=range_start;i<=range_end;i++){if(i==0)continue;print a[i],"spe"}  for(i=range_start;i<=range_end;i++){if(i==0)continue;print b[i],"uni"}  for(i=range_start;i<=range_end;i++){if(i==0)continue;print c[i],"other"} } ' ./enh_spe_H3K27ac_signal  ./enh_uni_H3K27ac_signal  ./enh_other_H3K27ac_signal >./enh_H3K27ac_signal

rm ./enh_spe_H3K4me1_signal  ./enh_uni_H3K4me1_signal  ./enh_other_H3K4me1_signal ./enh_spe_H3K4me3_signal  ./enh_uni_H3K4me3_signal  ./enh_other_H3K4me3_signal  ./enh_spe_H3K27ac_signal  ./enh_uni_H3K27ac_signal  ./enh_other_H3K27ac_signal 


###########################################
#生成wig格式的文件倒入igv即可。
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do
  #各个组织H3K4me1 spe的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=255,0,0\nfixedStep chrom=chr1 start=1 setp=1 span=1"} { if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if(NR>1&&$NF=="spe")print $tissue_id }'   enh_H3K4me1_signal >./histone_wig/H3K4me1/spe/$i.wig
  #各个组织H3K4me1 uni的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=255,0,0\nfixedStep chrom=chr1 start=1 setp=1 span=1"} {  if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="uni")print $tissue_id }'   enh_H3K4me1_signal >./histone_wig/H3K4me1/uni/$i.wig
  #各个组织H3K4me1 other的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=255,0,0\nfixedStep chrom=chr1 start=1 setp=1 span=1"} {  if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="other")print $tissue_id }'   enh_H3K4me1_signal >./histone_wig/H3K4me1/other/$i.wig

##########
  #各个组织H3K4me3 spe的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=0,0,255\nfixedStep chrom=chr1 start=1 setp=1 span=1"} {  if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="spe")print $tissue_id }'   enh_H3K4me3_signal >./histone_wig/H3K4me3/spe/$i.wig
  #各个组织H3K4me3 uni的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=0,0,255\nfixedStep chrom=chr1 start=1 setp=1 span=1"} {  if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="uni")print $tissue_id }'   enh_H3K4me3_signal >./histone_wig/H3K4me3/uni/$i.wig
  #各个组织H3K4me3 other的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=0,0,255\nfixedStep chrom=chr1 start=1 setp=1 span=1"} {  if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="other")print $tissue_id }'   enh_H3K4me3_signal >./histone_wig/H3K4me3/other/$i.wig
##########
  #各个组织H3K27ac spe的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=0,255,0\nfixedStep chrom=chr1 start=1 setp=1 span=1"} { if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="spe")print $tissue_id }'   enh_H3K27ac_signal >./histone_wig/H3K27ac/spe/$i.wig
  #各个组织H3K27ac uni的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=0,255,0\nfixedStep chrom=chr1 start=1 setp=1 span=1"} {  if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="uni")print $tissue_id }'   enh_H3K27ac_signal >./histone_wig/H3K27ac/uni/$i.wig
  #各个组织H3K27ac other的信号：
  awk -v tissue=$i 'BEGIN{OFS="";print "browser position chr1:1-4000\ntrack type=wiggle_0 name=\"",tissue,"\" color=0,255,0\nfixedStep chrom=chr1 start=1 setp=1 span=1"} {  if(NR==1){for(i=1;i<=NF;i++){if($i==tissue){tissue_id=i}}} if($12=="other")print $tissue_id }'   enh_H3K27ac_signal >./histone_wig/H3K27ac/other/$i.wig

echo $i
done













