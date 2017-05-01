#从shell读入工作目录
rm(list=ls())
now_path<-commandArgs()
now_path<-now_path[length(now_path)]
print(now_path)
setwd(now_path)

library(ggplot2)
library(reshape)


mean_signal_500<-as.data.frame(read.table("./enh_find_1/result.mean_signal",header=F,sep="\t",stringsAsFactors = F))
apply(mean_signal_500,2,as.numeric)
names(mean_signal_500)<-c("loc","H3K4me1","H3K4me3","H3K27ac")

# H3K4me1_sigsum<-as.numeric(read.table("./bg/H3K4me1_sigsum"))
# H3K4me3_sigsum<-as.numeric(read.table("./bg/H3K4me3_sigsum"))
# H3K27ac_sigsum<-as.numeric(read.table("./bg/H3K27ac_sigsum"))
# # # 
# mean_signal_500$H3K4me1<-mean_signal_500$H3K4me1*10^6/H3K4me1_sigsum
# mean_signal_500$H3K4me3<-mean_signal_500$H3K4me3*10^6/H3K4me3_sigsum
# mean_signal_500$H3K27ac<-mean_signal_500$H3K27ac*10^6/H3K27ac_sigsum

mean_signal_500_melt<-melt(mean_signal_500,id=c("loc"))
names(mean_signal_500_melt)[2]<-c("histone")
png_path<-paste(now_path,"/result.mean_signal.png",sep="")
png(png_path)
ggplot(data=mean_signal_500_melt,aes(x=loc,y=value,colour=histone))+geom_point()+scale_color_manual(values=c("blue","green","red"))
dev.off()
print(getwd())