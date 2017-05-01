#获得每个组织每个增强子上面的组蛋白修饰的信号平均值：
#定义函数来查看文件夹是否存在，不存在则建立：
function histone_dir_exists(){
    if [ ${1} != '' ]&&[ ! -d ./${1} ]; then
        echo ${1}文件夹及子文件夹不存在，已经建立
        mkdir -p ${1}/bed ${1}/signal
    fi
}
#使用函数创建：H3K4me1/H3K4me3/H3K27ac的文件夹及子文件夹：
cd /home/ding/all_cmd/script/enh_statistics
histone_dir_exists H3K4me1
histone_dir_exists H3K4me3
histone_dir_exists H3K27ac

###############

#定义函数：
function histone_signal(){
    if [  ${1} != '' ]; then
        cd /home/ding/all_cmd/script/enh_statistics
        tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)

        for i in ${tissue[@]}
        do 
            #先获得每个组织的enh的bed文件：
            awk -v tissue=$i 'BEGIN{OFS="\t"}{ if(NR>1&&$8==tissue){printf "%s\t%d\t%d\t%s\n",$1,$2,$3,$4} }' ./enh_all_count_spe_melt|sortBed -i - >./${1}/bed/$i"_enh.bed"
            echo $i"完成1"

            #通过bwtools获得每个增强子上的平均信号：
            bwtool summary ./${1}/bed/$i"_enh.bed" /media/ding/000B49000006264C/eRNA_project/histone/Adrenal/bg/${1}.bw -with-sum  -fill=0  /dev/stdout |awk -v tissue="$i" 'BEGIN{OFS="\t";FS="\t";print "ID","mean_signal","tissue"}{print $1":"$2"-"$3,$8,tissue}' - >./${1}/signal/$i"_enh_signal" 
            echo $i"完成2"

            ##添加增强子的名称到结果文件的最后一列：
            awk 'BEGIN{OFS="\t"}{ if(NR==FNR){a[$1":"$2"-"$3]=$4; } if(NR>FNR&&FNR==1){print $0,"enh_name"} if(NR>FNR&&FNR>1){print $0,a[$1]} }' ./${1}/bed/$i"_enh.bed" ./${1}/signal/$i"_enh_signal">./${1}/signal/a

            awk '{print $0}' ./${1}/signal/a >./${1}/signal/$i"_enh_signal"
            rm ./${1}/signal/a   
        done


        #将十个组织组蛋白信号合并：
        cd /home/ding/all_cmd/script/enh_statistics
        #添加列名并创建文件：
        awk 'BEGIN{OFS="\t";print "ID","mean_signal","tissue","enh_name"}' >./${1}/${1}.tsv
        #开始添加内容：
        tissue=(Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney)
        for i in ${tissue[@]}
            do 
            awk '{if(NR>1){print $0}}' ./${1}/signal/$i"_enh_signal" >>./${1}/${1}.tsv
            echo $i"完成4"
        done
    fi  
}

#获得H3K4me1/H3K4me3/H3K27ac的每个组织的每个增强子的信号值：
histone_signal H3K4me1
histone_signal H3K4me3
histone_signal H3K27ac


