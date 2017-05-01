#从shell读入工作目录
rm(list=ls())
library(ggplot2)
library(reshape)
library(stringr)
setwd("/home/ding/all_cmd/script")

spe_TF_enh_tissue<-read.table("./enh_statistics/a",header = F,stringsAsFactors = F)
#保留唯一的enh
spe_enh<-unique(spe_TF_enh_tissue$V9)
#创建新表格保留结果:
spe_TF_enh_tissue2<-unique(spe_TF_enh_tissue[,c("V9","V10","V14")])
names(spe_TF_enh_tissue2)<-c("enh","type","tissues")
for(i in spe_enh){
  spe_TF_enh_tissue2[spe_TF_enh_tissue2$enh==i,"TF"]<-paste(unique(spe_TF_enh_tissue[spe_TF_enh_tissue$V9==i,"V4"]),collapse = ",")
}

write.table(spe_TF_enh_tissue2,"./enh_statistics/spe_TF_enh_tissue_1",quote = F,col.names = F,row.names = F,sep="\t")






