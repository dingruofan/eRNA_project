#从shell读入工作目录
rm(list=ls())
library(ggplot2)
library(reshape)
library(stringr)
setwd("/home/ding/all_cmd/script")

spe_TF_enh_tissue<-read.table("/home/ding/all_cmd/script/enh_statistics/enh_TF/enh_TF_all",header = T,stringsAsFactors = F)
#保留唯一的enh
spe_enh<-unique(spe_TF_enh_tissue$enh_name)
#创建新表格保留结果:
spe_TF_enh_tissue2<-unique(spe_TF_enh_tissue[,c("enh_name","tissue")])
for(i in spe_enh){
  spe_TF_enh_tissue2[spe_TF_enh_tissue2$enh_name==i,"TF_name"]<-paste(unique(spe_TF_enh_tissue[spe_TF_enh_tissue$enh_name==i,"TF_name"]),collapse = ",")
}
spe_TF_enh_tissue2$type<-"-"
spe_TF_enh_tissue2<-spe_TF_enh_tissue2[,c("enh_name","type","TF_name","tissue")]

#再与TS enh取交集获得TS TF-enh的:
ts_enh<-read.table("/home/ding/all_cmd/script/enh_statistics/enh_spe_cluster.bed",header = F,stringsAsFactors = F)
names(ts_enh)<-c("chr","chr_start","chr_end","enh_name")
ts_enh<-ts_enh[,4,drop=F]

#合并:
spe_TF_enh_tissue3<-merge(spe_TF_enh_tissue2,ts_enh,by="enh_name")

write.table(spe_TF_enh_tissue3,"./enh_statistics/all_TF_enh_tissue_1",quote = F,col.names = F,row.names = F,sep="\t")






