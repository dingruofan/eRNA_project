rm(list=ls())
library(ggplot2)
library(corrplot)
#读入组蛋白修饰数据/DHS数据/GC含量数据
setwd("/home/ding/all_cmd/script/enh_statistics")

H3K4me1<-read.table("./H3K4me1/H3K4me1.tsv",header = T,stringsAsFactors = F)
H3K4me1<-H3K4me1[,c("mean_signal","tissue","enh_name")]
names(H3K4me1)[1]<-"H3K4me1"

H3K4me3<-read.table("./H3K4me3/H3K4me3.tsv",header = T,stringsAsFactors = F)
H3K4me3<-H3K4me3[,c("mean_signal","tissue","enh_name")]
names(H3K4me3)[1]<-"H3K4me3"

H3K27ac<-read.table("./H3K27ac/H3K27ac.tsv",header = T,stringsAsFactors = F)
H3K27ac<-H3K27ac[,c("mean_signal","tissue","enh_name")]
names(H3K27ac)[1]<-"H3K27ac"

GC<-read.table("./GC/GC_all.tsv",header = T,stringsAsFactors = F)
GC<-GC[GC$tissue!="random",]
names(GC)[2]<-"mean_signal"
GC<-GC[,c("mean_signal","tissue","enh_name")]
names(GC)[1]<-"GC"

DHS<-read.table("./DHS/DHS.tsv",header = T,stringsAsFactors = F)
DHS<-DHS[,c("mean_signal","tissue","enh_name")]
names(DHS)[1]<-"DHS"

#进行数据融合
a<-merge(H3K4me1,H3K27ac,by=c("tissue","enh_name"))
a<-merge(a,H3K4me3,by=c("tissue","enh_name"))
a<-merge(a,GC,by=c("tissue","enh_name"))
a<-merge(a,DHS,by=c("tissue","enh_name"))

#各个组织所有增强子相关性分析:
M<-cor(as.matrix(a[,c(3:7)]))
corrplot(M, method="circle")
corrplot(M, method="square")
corrplot(M, method="color")
corrplot(M, method="color",type="upper")
corrplot.mixed(M)
corrplot.mixed(M, lower="square", upper="circle")

#所有组织的相关性:
cor(as.matrix(a[,c(3:7)]))
#可以看出H3k27ac与GC含量相关性还可以0.3，H3K4me1和H3K4me3的相关性还可以0.42，其他的不行

#做每个组织的相关性分析:
tissues<-unique(H3K4me1$tissue)
for(i in tissues){
  cor_value<-cor.test(a[a$tissue==i,"GC"],a[a$tissue==i,"H3K27ac"])$estimate
  print(paste(i,cor_value,sep=":"))
}




