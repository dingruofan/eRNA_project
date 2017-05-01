#从shell读入工作目录
rm(list=ls())
args<-commandArgs()
now_path<-args[length(args)-1]
figure_path<-args[length(args)]
tissue<-args[length(args)-2]

# setwd("/media/ding/000B49000006264C/eRNA_project/histone/SkeletalMuscle")
library(ggplot2)
library(reshape)
library(Cairo)

mean_signal_500<-as.data.frame(read.table("./enh_find/result.mean_signal",header=F,sep="\t",stringsAsFactors = F))
apply(mean_signal_500,2,as.numeric)
names(mean_signal_500)<-c("loc","H3K4me1","H3K4me3","H3K27ac","Dnase I")

mean_signal_500_melt<-melt(mean_signal_500,id=c("loc"))
names(mean_signal_500_melt)[2]<-c("nor_signal")
png_path<-paste(figure_path,"_histone_DHS.png",sep="")
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 6, height = 6, units='in', dpi=600)

p<-ggplot(data=mean_signal_500_melt,aes(x=loc,y=value,colour=nor_signal))+
  geom_point()+
  scale_color_manual(values=c("blue","green","red","yellow"))+
  xlab("Distance to enhancer center (bp)")+
  ylab("Normalized mean signal")+
  annotate("text",x=min(mean_signal_500_melt$loc)+50,y=Inf,label=tissue,vjust=2.5,size=7,colour="#000080",hjust=0)+
  theme_bw()+
  theme(axis.text=element_text(size = rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        legend.text=element_text(size=rel(1.3)),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.margin=margin(0,0,0,0,"mm"),
        #legend.position=c(1,1),legend.justification=c(1,1),
        legend.position="top"
        )
plot(p)
dev.off()
