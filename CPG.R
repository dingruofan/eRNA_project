#从shell读入工作目录
rm(list=ls())
library(ggplot2)
library(Cairo)
setwd("/home/ding/all_cmd/script/enh_statistics")
#增强子GC:
CPG_dist<-read.table("CPG_dist",header = F,stringsAsFactors = F)
#与gene TSS比较：
CPG_gene_dist<-read.table("CPG_protein_gene2_dist",header = F,stringsAsFactors = F)



names(CPG_dist)<-c("enh","type","isspe","dist")
names(CPG_gene_dist)<-c("enh","type","isspe","dist")

CPG_dist$class<-"enh"
CPG_gene_dist$class<-"gene"

CPG_all_dist<-rbind(CPG_dist,CPG_gene_dist)

CPG_all_dist_table<-data.frame(t(table(CPG_all_dist$dist,CPG_all_dist$class)),stringsAsFactors = F)
names(CPG_all_dist_table)<-c("class","dist","Freq")
CPG_all_dist_table$dist<-as.integer(CPG_all_dist_table$dist)
CPG_all_dist_table[CPG_all_dist_table$Freq==0,"Freq"]<-NA
CPG_all_dist_table<-na.omit(CPG_all_dist_table)

ggplot(CPG_all_dist_table,aes(x=dist,y=Freq))+geom_point()+xlim(-2000,2000)+ylim(0,100)

binwidth_2<-diff(range(CPG_all_dist$dist))/5000000
ggplot(CPG_all_dist,aes(dist))+
  geom_histogram(binwidth=binwidth_2,aes(fill=class))+xlim(-2000,2000)+ylim(0,50)
  geom_density(aes(fill=class))+xlim(-2000,2000)

# ggplot(CPG_dist,aes(x=dist,y=..density..))+geom_density(aes(colour=type))+xlim(-2000,2000)
# 
# ggplot(CPG_dist,aes(x=dist,y=..density..))+geom_density(aes(colour=isspe))+xlim(-2000,2000)
  
png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/CPG_2.png"
# 使用Cairo生成高清图片：par('din')获得当前显示图的尺寸：
CairoPNG(png_path, width = 6.5, height = 6, units='in', dpi=600)

binwidth_1<-diff(range(CPG_dist$dist))/80000

ggplot(CPG_dist,aes(x=dist,y=..density..))+
  geom_histogram(binwidth = binwidth_1,alpha=0.3)+
  geom_density(alpha=0.8)+
  xlim(-2000,2000)+
  xlab("distance to enhancer(bp)")+
  theme_bw()+
  theme(axis.text=element_text(size=rel(1.0)),
        axis.text.x=element_text(vjust=1),
        axis.title=element_text(size=rel(1.1)),
        legend.text=element_text(size=rel(1.0)),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        # legend.margin=unit(0,"mm"),
        #legend.position=c(1,1),legend.justification=c(1,1),
        legend.position="top"
  )

dev.off()

CPG_dist_overlap<-CPG_dist[CPG_dist$dist==0,]
table(CPG_dist_overlap$isspe)
#共有982, other 264   spe 688    uni30 

table(cut(CPG_dist$dist,breaks=c(-10^(6:4),-2000,-1,0,2000,10^(4:6)),labels = paste(c("-1000~-100","-100~-10","-10~-2","-2~-0.001","-0.001~0","0~2","2~10","10~100","100~1000"),"kb",sep="")))
table(cut(CPG_all_dist$dist,breaks=c(-10^(6:4),-2000,-1,0,2000,10^(4:6)),labels = paste(c("-1000~-100","-100~-10","-10~-2","-2~-0.001","-0.001~0","0~2","2~10","10~100","100~1000"),"kb",sep="")))



########################
#做CPG的bigwig在增强子上的覆盖:
rm(list=ls())
library(ggplot2)
library(Cairo)
setwd("/home/ding/all_cmd/script/enh_statistics")
#增强子GC:
CPG_dist<-read.table("CPG_enh_dist0.mean_signal",header = F,stringsAsFactors = F)
names(CPG_dist)<-c("location","CpG_signal")

#注意这里:！!!！！!！!！!！
##不需要计算平均信号，所以每个位置都乘982
CPG_dist$CpG_signal_2<-CPG_dist$CpG_signal*640


# # 平滑：
# i<--2000
# window_size=10
# mean_signal_500_smooth<-data.frame(loc=seq(-2000,2000,1))
# while(i<2000){
#   mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"random"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"random"])
#   mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"2D-eRNAs Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"2D-eRNAs Enh"])
#   mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"1D-eRNAs Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"1D-eRNAs Enh"])
#   mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"Gene TSS"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"Gene TSS"])
#   mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"no-eRNAs Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"no-eRNAs Enh"])
#   
#   i=i+window_size
# }
# #只能用一个
# mean_signal_500_melt<-melt(mean_signal_500_smooth,id=c("loc"))
# #只能用一个
# mean_signal_500_melt<-melt(mean_signal_500,id=c("loc"))


png_path="/media/ding/000B49000006264C/eRNA_project/figure/CpG_signal.png"
CairoPNG(png_path, width = 6, height = 5.5, units='in', dpi=600)

# png_path="/media/ding/000B49000006264C/eRNA_project/figure/GC_isbid.png"
# CairoPNG(png_path, width = 6, height = 5.5, units='in', dpi=600)

ggplot(data=CPG_dist,aes(x=location,y=CpG_signal_2))+
  geom_line(size=0.9)+
  # geom_smooth(se=FALSE,aes(fill=type))+
  # scale_color_manual(values=c("blue","green","red","black"))+
  xlab("Distance to enhancer center (bp)")+
  ylab("Signal of CpG island ")+
  # geom_vline(xintercept = -100)+
  scale_x_continuous(breaks=seq(-2000,2000,1000),labels=paste(seq(-2000,2000,1000),seq=""))+
  
  theme_bw()+
  theme(
    axis.text=element_text(size=rel(1.1)),
        axis.title=element_text(size=rel(1.1)),
        legend.text=element_text(size=rel(1.1)),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.margin=margin(0,0,0,0,"pt"),
        #legend.position=c(1,1),legend.justification=c(1,1),
        legend.position="top",
        legend.box.spacing = unit(0,"mm")
  )

dev.off()



