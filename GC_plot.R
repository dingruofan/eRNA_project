#???shell读???工?????????
# rm(list=ls())
# now_path<-commandArgs()
# now_path<-now_path[length(now_path)]
# print(now_path)
# setwd(now_path)

rm(list=ls())
setwd("/home/ding/all_cmd/script")

library(ggplot2)
library(reshape)
library(Cairo)
options(bitmapType = "cairo")

mean_signal_500<-as.data.frame(read.table("./enh_statistics/GC_result.mean_signal",header=F,sep="\t",stringsAsFactors = F))
apply(mean_signal_500,2,as.numeric)
names(mean_signal_500)<-c("loc","random","2D-eRNAs Enh","1D-eRNAs Enh","Gene TSS","no-eRNAs Enh")

# 平滑：
i<--2000
window_size=10
mean_signal_500_smooth<-data.frame(loc=seq(-2000,2000,1))
while(i<2000){
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
mean_signal_500_melt$type<-factor(mean_signal_500_melt$type,levels=c("2D-eRNAs Enh","1D-eRNAs Enh","no-eRNAs Enh","Gene TSS","random"))
  
  #画图

  png_path="/media/ding/000B49000006264C/eRNA_project/figure/GC_isbid.png"
  CairoPNG(png_path, width = 6.2, height = 5.5, units='in', dpi=600)
  
  # png_path="/media/ding/000B49000006264C/eRNA_project/figure/GC_isbid.png"
  # CairoPNG(png_path, width = 6, height = 5.5, units='in', dpi=600)
  
  ggplot(data=mean_signal_500_melt,aes(x=loc,y=value,colour=type))+
    geom_line(size=0.6)+
    # geom_smooth(se=FALSE,aes(fill=type))+
    # scale_color_manual(values=c("blue","green","red","black"))+
    xlab("Distance to enhancer center/ TSS of gene (bp)")+
    ylab("GC %")+
    # geom_vline(xintercept = -100)+
    scale_colour_manual(limits=c("2D-eRNAs Enh","1D-eRNAs Enh","no-eRNAs Enh","Gene TSS","random"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"]),"Gene TSS","random"),values = c("red","gold","green","blue","black"))+
    scale_x_continuous(breaks=seq(-2000,2000,1000),labels=paste(seq(-2000,2000,1000),seq=""))+
  
    theme_bw()+
    theme(
      axis.text=element_text(size = rel(1.3)),
      axis.title=element_text(size=rel(1.3)),
      # axis.title.x=element_blank(),
      legend.text=element_text(size=rel(1.2)),
      legend.title=element_blank(),
      legend.background=element_blank(),
      legend.key = element_blank(),
      legend.margin=margin(0,0,0,0,"mm"),
      legend.box.spacing = unit(1,"mm"),
      #legend.position=c(1,1),legend.justification=c(1,1),
      legend.position="top"
    )
    
dev.off()

#每个组织增强子的GC含量boxplot:

rm(list=ls())
setwd("/home/ding/all_cmd/script")

library(ggplot2)
library(reshape)
library(Cairo)
options(bitmapType = "cairo")

GC_all<-as.data.frame(read.table("./enh_statistics/./GC/GC_all.tsv",header=T,sep="\t",stringsAsFactors = F))
GC_all$GC_percent<-as.numeric(GC_all$GC_percent)

#标注组织和random不同的颜色:
GC_all[GC_all$tissue=="random","colours"]<-"grey"
GC_all[GC_all$tissue!="random","colours"]<-"blue"
# cols<-c("1"="red","0"="green")

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/GC_boxplot.png"
CairoPNG(png_path, width = 5.5, height =5.16 , units='in', dpi=600)

ggplot(GC_all,aes(x=reorder(tissue,GC_percent,median),y=GC_percent))+
  geom_boxplot(aes(fill=colours),show.legend = T,outlier.size = 0,outlier.stroke = 0)+
  # ylim(20,80)+
  xlab("tissue")+
  ylab("GC content(%)")+
  scale_fill_manual(values = c("#FF5151","grey"))+
  scale_y_continuous(breaks=seq(0,80,20))+
  theme_bw()+
  theme(axis.text=element_text(size=rel(1.3)),
        axis.text.x=element_text(angle = 30,vjust=1,hjust=1),
        legend.text=element_text(size=rel(1.3)),
        axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.margin=margin(0,0,0,0,"mm"),
        #legend.position=c(1,1),legend.justification=c(1,1),
        legend.position="top",
        legend.box.spacing = unit(0,"mm")
  )+
  guides(fill=F)

dev.off()

t.test(GC_all[GC_all$tissue=="Breast",]$GC_percent,GC_all[GC_all$tissue=="random",]$GC_percent,alternative = "greater")


##########
#临时计算p值，无影响
rm(list=ls())
setwd("/home/ding/all_cmd/script/enh_statistics/GC/GC_percent")
tissues_GC<-read.table("all_enh_GC",header = T,stringsAsFactors = F)
random_GC<-read.table("random_GC",header = T,stringsAsFactors = F)

t.test(tissues_GC$GC_percent,random_GC$GC_percent,alternative = "greater")$p.value

#加入每个增强子的注释：
rm(list=ls())
setwd("/home/ding/all_cmd/script/enh_statistics/GC/GC_percent")
GC_all<-read.table("../GC_all.tsv",header = T,stringsAsFactors = F)
GC_all_enh<-GC_all[GC_all$tissue!="random",]
GC_all_random<-GC_all[GC_all$tissue=="random",]
GC_all_enh$tissue<-NULL
GC_all_enh<-GC_all_enh[!duplicated(GC_all_enh),]
#读入增强子数据:
enh_all<-read.table("/home/ding/all_cmd/script/enh_statistics/enh_all",header = T,stringsAsFactors = F)
names(enh_all)[4]<-"enh_name"


all_GC_type<-merge(GC_all_enh,enh_all,by="enh_name")

t.test(all_GC_type$GC_percent,GC_all_random$GC_percent,alternative = "greater")$p.value
#2D>1D:
t.test(all_GC_type[all_GC_type$type=="bid",]$GC_percent,all_GC_type[all_GC_type$type=="unbid",]$GC_percent,alternative = "greater")$p.value
#1D>NO
t.test(all_GC_type[all_GC_type$type=="unbid",]$GC_percent,all_GC_type[all_GC_type$type=="no_eRNA",]$GC_percent,alternative = "greater")$p.value


#UE>Oth:
t.test(all_GC_type[all_GC_type$is_spe=="uni",]$GC_percent,all_GC_type[all_GC_type$is_spe=="other",]$GC_percent,alternative = "greater")$p.value
#Oth>TS
t.test(all_GC_type[all_GC_type$is_spe=="other",]$GC_percent,all_GC_type[all_GC_type$is_spe=="spe",]$GC_percent,alternative = "greater")$p.value


#eRNA >no_eRNA
t.test(all_GC_type[all_GC_type$type=="bid"|all_GC_type$type=="unbid",]$GC_percent,all_GC_type[all_GC_type$type=="no_eRNA",]$GC_percent,alternative = "greater")$p.value





