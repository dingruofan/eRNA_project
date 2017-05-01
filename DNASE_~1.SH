#!/bin/bash
#先把10个组织的DHS进行merge然后再分配给各个组织：

cd /media/ding/000B49000006264C/eRNA_project/DNase
bedops --merge ./Adrenal/Dnase_ENCODE.bed ./Brain/Dnase_ENCODE.bed ./Breast/Dnase_ENCODE.bed ./Heart/Dnase_ENCODE.bed ./Liver/Dnase_ENCODE.bed ./Lung/Dnase_ENCODE.bed ./Ovary/Dnase_ENCODE.bed ./Placenta/Dnase_ENCODE.bed ./SkeletalMuscle/Dnase_ENCODE.bed ./Kidney/Dnase_ENCODE.bed  >/home/ding/all_cmd/script/enh_statistics/Dnase_cluster.bed

#分配给各个组织：
cd /media/ding/000B49000006264C/eRNA_project/DNase
bedops --element-of 1  /home/ding/all_cmd/script/enh_statistics/Dnase_cluster.bed ./Adrenal/Dnase_ENCODE.bed >./Adrenal/Dnase_ENCODE_cluster.bed
bedops --element-of 1  /home/ding/all_cmd/script/enh_statistics/Dnase_cluster.bed ./Brain/Dnase_ENCODE.bed >./Brain/Dnase_ENCODE_cluster.bed
bedops --element-of 1  /home/ding/all_cmd/script/enh_statistics/Dnase_cluster.bed ./Breast/Dnase_ENCODE.bed >./Breast/Dnase_ENCODE_cluster.bed
bedops --element-of 1  /home/ding/all_cmd/script/enh_statistics/Dnase_cluster.bed ./Heart/Dnase_ENCODE.bed >./Heart/Dnase_ENCODE_cluster.bed
bedops --element-of 1  /home/ding/all_cmd/script/enh_statistics/Dnase_cluster.bed ./Liver/Dnase_ENCODE.bed >./Liver/Dnase_ENCODE_cluster.bed
bedops --element-of 1  /home/ding/all_cmd/script/enh_statistics/Dnase_cluster.bed ./Lung/Dnase_ENCODE.bed >./Lung/Dnase_ENCODE_cluster.bed
bedops --element-of 1  /home/ding/all_cmd/script/enh_statistics/Dnase_cluster.bed ./Ovary/Dnase_ENCODE.bed >./Ovary/Dnase_ENCODE_cluster.bed
bedops --element-of 1  /home/ding/all_cmd/script/enh_statistics/Dnase_cluster.bed ./Placenta/Dnase_ENCODE.bed >./Placenta/Dnase_ENCODE_cluster.bed
bedops --element-of 1  /home/ding/all_cmd/script/enh_statistics/Dnase_cluster.bed ./SkeletalMuscle/Dnase_ENCODE.bed >./SkeletalMuscle/Dnase_ENCODE_cluster.bed
bedops --element-of 1  /home/ding/all_cmd/script/enh_statistics/Dnase_cluster.bed ./Kidney/Dnase_ENCODE.bed >./Kidney/Dnase_ENCODE_cluster.bed

echo "*********************************************"
echo -e "\tDHS 聚类 OK!"
echo "*********************************************"


