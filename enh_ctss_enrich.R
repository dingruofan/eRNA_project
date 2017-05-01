# #从shell读入工作目录
# rm(list=ls())
# now_path<-commandArgs()
# now_path<-now_path[length(now_path)]
# print(now_path)
# setwd(now_path)

rm(list=ls())
setwd("/media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/")

library(ggplot2)
library(reshape2)
library(Cairo)
#17537

mean_signal_500<-as.data.frame(read.table("./result.mean_signal",header=F,sep="\t",stringsAsFactors = F))
apply(mean_signal_500,2,as.numeric)
names(mean_signal_500)<-c("loc","pos","neg")
# mean_signal_500$pos<-mean_signal_500$pos/sum(mean_signal_500$pos)
# mean_signal_500$neg<-mean_signal_500$neg/sum(mean_signal_500$neg)
#用窗口平滑:
# i=-5000
# mean_signal_win<-data.frame(loc=seq(-5000,5000,1),pos="",neg="",stringsAsFactors = F)
# while(i<=5000){
#   mean_signal_win[which(mean_signal_win$loc %in% seq(i,i+50,1)),"pos"]<-mean(mean_signal_500[which(mean_signal_500$loc %in% seq(i,i+50,1)),"pos"])
#   mean_signal_win[which(mean_signal_win$loc %in% seq(i,i+50,1)),"neg"]<-mean(mean_signal_500[which(mean_signal_500$loc %in% seq(i,i+50,1)),"neg"])
#   i=i+50
# }
mean_signal_win<-mean_signal_500

mean_signal_500_melt<-melt(mean_signal_win,id=c("loc"))
mean_signal_500_melt$value<-as.numeric(mean_signal_500_melt$value)
names(mean_signal_500_melt)[2]<-c("strand")
ggplot(data=mean_signal_500_melt,aes(x=loc,y=value,colour=strand))+geom_smooth(position = "identity",se=FALSE,n=10000,span=0)+scale_color_manual(values=c("red","blue"))+theme(axis.text=element_text(size=20),legend.text=element_text(size=20))+ylim(0,0.02)+xlim(-500,500)
ggplot(data=mean_signal_500_melt,aes(x=loc,y=value,colour=strand))+geom_smooth(position = "identity",se=FALSE,n=10000,span=0)+scale_color_manual(values=c("red","blue"))+theme(axis.text=element_text(size=20),legend.text=element_text(size=20))+ylim(0,0.1)+xlim(-5000,5000)

mean_signal_500_melt<-melt(mean_signal_win,id=c("loc"))
mean_signal_500_melt$value<-as.numeric(mean_signal_500_melt$value)
names(mean_signal_500_melt)[2]<-c("strand")



#500bp：
png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_ctss_enrich_500.png"
CairoPNG(png_path, width = par('din')[1], height = par('din')[2], units='in', dpi=600)

ggplot(data=mean_signal_500_melt,aes(x=loc,y=value,colour=strand))+
  geom_smooth(position = "identity",se=FALSE,n=10000,span=0)+
  scale_color_manual(values=c("red","blue"))+
  ylim(0,0.02)+xlim(-500,500)+
  xlab("distance to enhancer center (bp)")+
  ylab("density")+
  scale_colour_discrete(limits=c("pos","neg"),labels=c("CAGE +  ","CAGE -  "))+
  theme_bw()+
  theme(
        legend.text=element_text(size=rel(1.0)),
        axis.text=element_text(size = rel(1.0)),
        axis.title=element_text(size=rel(1.0)),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.margin=unit(0.1,"mm"),
        #legend.position=c(1,1),legend.justification=c(1,1),
        legend.position="top"
        )
  
dev.off()


png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_ctss_enrich_5000.png"
CairoPNG(png_path, width = 6.7, height = 6.4, units='in', dpi=600)

ggplot(data=mean_signal_500_melt,aes(x=loc,y=value,colour=strand))+
  geom_smooth(position = "identity",se=FALSE,n=10000,span=0.9)+
  geom_point(alpha=0.1,size=0.5)+
  scale_color_manual(values=c("red","blue"))+
  ylim(0,0.1)+xlim(-5000,5000)+
  xlab("Distance to enhancer center (kb)")+
  ylab("Mean signal")+
  scale_colour_discrete(limits=c("pos","neg"),labels=c("CTSS +  ","CTSS -  "))+
  scale_x_continuous(breaks=seq(-5000,5000,1000),labels=seq(-5,5,1),limits = c(-5000,5000))+
  theme_bw()+
  theme(
        legend.text=element_text(size=rel(1.0)),
        axis.text=element_text(size = rel(1.0)),
        axis.title=element_text(size=rel(1.0)),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.margin=margin(0,0,0,0,"mm"),
        #legend.position=c(1,1),legend.justification=c(1,1),
        legend.position="top",
        legend.box.spacing = unit(0.1,"mm"),
        # legend.key.size = unit(50,"mm")
  )

dev.off()


############
#所有增强子与gene的距离：

rm(list=ls())
setwd("/home/ding/all_cmd/script/enh_statistics")

library(ggplot2)
library(reshape2)
library(Cairo)

enh_dist<-read.table("./enh_allgene_dist",header = F,stringsAsFactors = F)
enh_dist_1<-enh_dist[,c(4,5,6,10,11,12,14)]
names(enh_dist_1)<-c("strand","gene_id","gene_type","enh","type","is_spe","dist")
enh_dist_1[enh_dist_1$gene_type!="protein_coding","gene_type"]<-"non-coding"

enh_dist_1<-enh_dist_1[enh_dist_1$dist<5000&enh_dist_1$dist>-5000,]
enh_dist_1_table<-as.data.frame(table(enh_dist_1$dist,enh_dist_1$gene_type))
names(enh_dist_1_table)<-c("location","gene_type","Freq")

ggplot(enh_dist_1_table,aes(x=location,y=Freq))+geom_smooth(,aes(fill=gene_type))+ylim(0,10)


enh_dist_2<-enh_dist_1[enh_dist_1$dist<1000&enh_dist_1$dist>-1000&enh_dist_1$gene_type=="protein_coding",]

# 
# #######################
# rm(list=ls())
# setwd("/media/ding/000B49000006264C/eRNA_project/histone/liver/")
# 
# library(ggplot2)
# library(reshape2)
# 
# mean_signal_500<-as.data.frame(read.table("./result.mean_signal",header=F,sep="\t",stringsAsFactors = F))
# apply(mean_signal_500,2,as.numeric)
# names(mean_signal_500)<-c("loc","pos","neg")
# mean_signal_win<-mean_signal_500
# 
# mean_signal_500_melt<-melt(mean_signal_win,id=c("loc"))
# mean_signal_500_melt$value<-as.numeric(mean_signal_500_melt$value)
# names(mean_signal_500_melt)[2]<-c("strand")
# ggplot(data=mean_signal_500_melt,aes(x=loc,y=value,colour=strand))+geom_smooth(position = "identity",se=FALSE,n=10000,span=0)+scale_color_manual(values=c("red","blue"))+theme(axis.text=element_text(size=20),legend.text=element_text(size=20))+ylim(0,0.00)+xlim(-500,500)
# ggplot(data=mean_signal_500_melt,aes(x=loc,y=value,colour=strand))+geom_smooth(position = "identity",se=FALSE,n=10000,span=0)+scale_color_manual(values=c("red","blue"))+theme(axis.text=element_text(size=20),legend.text=element_text(size=20))+
#  xlim(-5000,5000)
# # 
# ggplot(data=mean_signal_500_melt,aes(x=loc,y=value,colour=strand))+geom_point(position = "identity")



# 
# 
# ################
# rm(list=ls())
# library(ggplot2)
# library(reshape2)
# 
# now_path<-commandArgs()
# now_path<-now_path[length(now_path)]
# print(now_path)
# setwd(now_path)
# 
# # setwd("/media/ding/000B49000006264C/eRNA_project/histone/Brain")
# 
# pos_dist<-read.table("./pos.dist",header = F,stringsAsFactors = F)
# neg_dist<-read.table("./neg.dist",header = F,stringsAsFactors = F)
# names(pos_dist)<-c("dist","Freq")
# names(neg_dist)<-c("dist","Freq")
# 
# all_dist<-data.frame(dist=seq(-500,500,1))
# 
# pos_dist<-merge(all_dist,pos_dist,all.x = T)
# pos_dist[is.na(pos_dist)]<-0
# 
# #用50的窗口平滑
# i=-500
# pos_dist_windows<-data.frame(dist=seq(-500,500,1),Freq="",stringsAsFactors = F)
# while(i<=500){
#   pos_dist_windows[pos_dist_windows$dist==i,"Freq"]<-mean(pos_dist[which(pos_dist$dist %in% seq(i,i+1,1)),"Freq"])
#   i=i+1
# }
# pos_dist_windows$Freq<-as.numeric(pos_dist_windows$Freq)
# pos_dist_windows$Freq<-pos_dist_windows$Freq/sum(pos_dist_windows$Freq)
# pos_dist_windows$type<-"pos"
# 
# #######neg
# 
# neg_dist<-merge(all_dist,neg_dist,all.x = T)
# neg_dist[is.na(neg_dist)]<-0
# 
# #用50的窗口平滑
# i=-500
# neg_dist_windows<-data.frame(dist=seq(-500,500,1),Freq="",stringsAsFactors = F)
# while(i<=500){
#   neg_dist_windows[neg_dist_windows$dist==i,"Freq"]<-mean(neg_dist[which(neg_dist$dist %in% seq(i,i+1,1)),"Freq"])
#   i=i+1
# }
# neg_dist_windows$Freq<-as.numeric(neg_dist_windows$Freq)
# neg_dist_windows$Freq<-neg_dist_windows$Freq/sum(neg_dist_windows$Freq)
# neg_dist_windows$type<-"neg"
# 
# all_dist<-rbind(pos_dist_windows,neg_dist_windows)
# 
# # ggplot(all_dist,aes(x=dist,y=Freq))+geom_smooth(aes(color=type),se=FALSE)+ylim(0,0.2)
# 
# png_path<-paste(now_path,"/result.mean_signal.png",sep="")
# png(png_path)
# ggplot(all_dist,aes(x=dist,y=Freq))+geom_smooth(aes(color=type),se=FALSE)+geom_point(aes(color=type))
# dev.off()












