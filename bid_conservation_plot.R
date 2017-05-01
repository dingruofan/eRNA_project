#从shell读入工作目录
rm(list=ls())
now_path<-commandArgs()
now_path<-now_path[length(now_path)]
print(now_path)
setwd(now_path)


library(ggplot2)
library(reshape)


mean_signal_500<-as.data.frame(read.table("./enh_find/result.mean_signal",header=F,sep="\t",stringsAsFactors = F))
apply(mean_signal_500,2,as.numeric)
names(mean_signal_500)<-c("loc","phastCons")


png_path<-paste(now_path,"/result.mean_signal.png",sep="")
png(png_path)
ggplot(data=mean_signal_500,aes(x=loc,y=phastCons))+geom_point()
dev.off()
