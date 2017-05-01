#获得每个组织每个增强子中心100bp的DHS信号：

#先获得每个组织的enh的bed文件：
#只针对中心100bp的序列：
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  awk -v tissue=$i 'BEGIN{OFS="\t"}{ if(NR>1&&$8==tissue){printf "%s\t%d\t%d\t%s\n",$1,($2+$3)/2-100,($2+$3)/2+100,$4} }' enh_all_count_spe_melt|sortBed -i - >./DHS/bed/$i"_enh.bed"
echo $i"完成1"
done



#获得DHS信号，将NA替换为0 -fill=0
cd /home/ding/all_cmd/script/enh_statistics
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
  bwtool summary ./DHS/bed/$i"_enh.bed"  /media/ding/000B49000006264C/eRNA_project/DNase/$i/new_DHS_TPM.bw -with-sum /dev/stdout -fill=0 |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$NF}' -| sort  -k 4 -n -r - |awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1"\t"$2"\t"$3]=$4} if(NR>FNR){print a[$1"\t"$2"\t"$3],$4} }' ./DHS/bed/$i"_enh.bed" -  >./DHS/signal/$i"_enh_DHS.bed"

echo $i"完成2"
done

#合并所有组织的DHS信号到一个文件：
cd /home/ding/all_cmd/script/enh_statistics
#添加列名并创建文件：
awk 'BEGIN{OFS="\t";print "ID","mean_signal","tissue","enh_name"}' >./DHS/DHS.tsv
#开始添加内容：
tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
for i in ${tissue[@]}
do 
    awk -v tissue=$i '{print NR,$2,tissue,$1}' ./DHS/signal/$i"_enh_DHS.bed" >>./DHS/DHS.tsv
    echo $i"完成4"
done
