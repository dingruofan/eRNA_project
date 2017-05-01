
rm(list=ls())
setwd("/home/ding/all_cmd/script")

library(ggplot2)
library(reshape)
library(Cairo)

mean_signal_500<-as.data.frame(read.table("./enh_statistics/GC_vista_result.mean_signal",header=F,sep="\t",stringsAsFactors = F))
apply(mean_signal_500,2,as.numeric)
names(mean_signal_500)<-c("loc","vista Enh","random","2D-eRNAs Enh","1D-eRNAs Enh","Gene TSS","no-eRNAs Enh")

# 平滑：
i<--2000
window_size=10
mean_signal_500_smooth<-data.frame(loc=seq(-2000,2000,1))
while(i<2000){
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"vista Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"vista Enh"])
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"random"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"random"])
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"2D-eRNAs Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"2D-eRNAs Enh"])
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"1D-eRNAs Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"1D-eRNAs Enh"])
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"Gene TSS"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"Gene TSS"])
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"no-eRNAs Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"no-eRNAs Enh"])
  
  i=i+window_size
}
#只能用一个
mean_signal_500_melt<-melt(mean_signal_500_smooth,id=c("loc"))
#只能用一个
# mean_signal_500_melt<-melt(mean_signal_500,id=c("loc"))

names(mean_signal_500_melt)[2]<-c("type")

#排序:
mean_signal_500_melt$type<-factor(mean_signal_500_melt$type,levels=c("vista Enh","2D-eRNAs Enh","1D-eRNAs Enh","no-eRNAs Enh","Gene TSS","random"))

#画图
png_path="/media/ding/000B49000006264C/eRNA_project/figure/GC_vista_isbid.png"
CairoPNG(png_path, width = 6, height = 5.5, units='in', dpi=600)

ggplot(data=mean_signal_500_melt,aes(x=loc,y=value,colour=type))+
  geom_line(size=0.6)+
  # geom_smooth(se=FALSE,aes(fill=type))+
  # scale_color_manual(values=c("blue","green","red","black"))+
  # xlab("distance to enhancer center (bp)")+
  ylab("GC enrich")+
  scale_colour_manual(limits=c("vista Enh","2D-eRNAs Enh","1D-eRNAs Enh","no-eRNAs Enh","Gene TSS","random"),labels=c("vista Enh",expression(En["2D-eRNA"]),expression(En["1D-eRNA"]),expression(En["no-eRNA"]),"gene TSS","random"),values = c("purple","red","gold","green","blue","black"))+
  scale_x_continuous(breaks=seq(-2000,2000,1000),labels=paste(seq(-2000,2000,1000),"bp",seq=""))+

  theme_bw()+
  theme(axis.text=element_text(size=rel(1.1)),
        legend.text=element_text(size=rel(1.1)),
        axis.text=element_text(size = rel(1.1)),
        axis.title=element_text(size=rel(1.1)),
        axis.title.x=element_blank(),
        legend.text=element_text(size=rel(1.1)),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.margin=unit(0,"mm"),
        #legend.position=c(1,1),legend.justification=c(1,1),
        legend.position="top"
  )
  
dev.off()



