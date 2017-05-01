#???shell读???工?????????
# rm(list=ls())
# now_path<-commandArgs()
# now_path<-now_path[length(now_path)]
# print(now_path)
# setwd(now_path)
# # 
rm(list=ls())
setwd("/home/ding/all_cmd/script/")


library(ggplot2)
library(reshape)
library(Cairo)

# setwd("/home/ding/all_cmd/script")
mean_signal_500<-as.data.frame(read.table("./enh_statistics/conservation2_result.mean_signal",header=F,sep="\t",stringsAsFactors = F))
apply(mean_signal_500,2,as.numeric)
names(mean_signal_500)<-c("loc","random","specific Enh","ubiquitously Enh","gene TSS","other Enh")

# 平滑：
i<--2000
window_size=10
mean_signal_500_smooth<-data.frame(loc=seq(-2000,2000,1))
while(i<2000){
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"random"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"random"])
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"specific Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"specific Enh"])
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"ubiquitously Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"ubiquitously Enh"])
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"gene TSS"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"gene TSS"])
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"other Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"other Enh"])
  
  i=i+window_size
}
#只能用一个
mean_signal_500_melt<-melt(mean_signal_500_smooth,id=c("loc"))
#只能用一个
# mean_signal_500_melt<-melt(mean_signal_500,id=c("loc"))
names(mean_signal_500_melt)[2]<-c("type")

#排序：
mean_signal_500_melt$type<-factor(mean_signal_500_melt$type,levels=c("specific Enh","ubiquitously Enh","other Enh","gene TSS","random"))

#画图:
png_path="/media/ding/000B49000006264C/eRNA_project/figure/conservation_isspe.png"
CairoPNG(png_path, width = 6.2, height = 6, units='in', dpi=600)

ggplot(data=mean_signal_500_melt,aes(x=loc,y=value,colour=type))+
  # geom_smooth()+
  geom_line(size=0.6,alpha=1)+
  # geom_point()+
  xlab("Distance to enhancer center/ TSS of gene (bp)")+
  ylab("PhastCons Score")+
  xlim(c(-1000,1000))+
  scale_colour_manual(limits=c("specific Enh","ubiquitously Enh","other Enh","gene TSS","random"),labels=c(expression(Enh["TS"]),expression(Enh["UE"]),expression(Enh["Oth"]),"Gene TSS","random"),values = c("red","gold","green","blue","black"))+
  theme_bw()+
  theme(
    axis.text=element_text(size = rel(1.3)),
    axis.title=element_text(size=rel(1.3)),
    # axis.title.x=element_blank(),
    legend.text=element_text(size=rel(1.3)),
    legend.title=element_blank(),
    legend.background=element_blank(),
    legend.key = element_blank(),
    legend.margin=margin(0,0,0,0,"mm"),
    legend.box.spacing = unit(1,"mm"),
    #legend.position=c(1,1),legend.justification=c(1,1),
    legend.position="top"
  )

dev.off()
