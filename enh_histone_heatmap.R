#从shell读入工作目录
rm(list=ls())
library(ggplot2)
library(reshape2)
setwd("/home/ding/all_cmd/script/enh_statistics")
#读入每种组蛋白的每个组织的平均信号值:
H3K4me1<-read.table("./enh_H3K4me1_signal",header = T,stringsAsFactors = F)
H3K4me1_melt<-melt(H3K4me1,id.vars = c("location","type"),variable.name = "tissue",value.name = "signal")
H3K4me3<-read.table("./enh_H3K4me3_signal",header = T,stringsAsFactors = F)
H3K27ac<-read.table("./enh_H3K27ac_signal",header = T,stringsAsFactors = F)

ggplot(H3K4me1_melt,aes(x=location,y=tissue))+ stat_density(aes(fill = type), geom = "raster", position = "identity")
heatmap3(as.matrix(H3K4me1[1:4000,2:11]))


ggplot(H3K4me1[,c("location","Adrenal","type")],aes(x=location,y=Adrenal))+geom_
