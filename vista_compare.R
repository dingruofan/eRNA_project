# rm(list=ls())
# library(ggplot2)
# now_path<-commandArgs()
# now_path<-now_path[length(now_path)]
# print(now_path)
# setwd(now_path)

setwd("/home/ding/all_cmd/script")

#统计vista enhancer的长度：
vista_enh<-read.table("/media/ding/000B49000006264C/eRNA_project/vista_enhancer/midbrain_mesencephalon/vista_enh_outgene.bed",header = F,stringsAsFactors = F)
vista_enh$length<-vista_enh$V3-vista_enh$V2
mean(vista_enh$length)

dist<-read.table("./enh_statistics/dist_to_plot",stringsAsFactors = F)
names(dist)<-c("enh_type","dist")
dist$dist<-abs(dist$dist)

#对距离进行分割:
dist_table<-data.frame(t(table(dist$dist)))[,c(2,3)]
names(dist_table)<-c("dist","Freq")

# dist_table_cut<-data.frame(t(table(cut(dist$dist,breaks=c(-1,1,10^3,10^4,10^5,10^6),labels = c("has Overlap","1-1kb","1kb-10kb","10kb-100kb",">100kb")))))[,c(2,3)]
dist_table_cut<-data.frame(t(table(cut(dist$dist,breaks=c(-1,1,10^5,10^6),labels = c("has Overlap","1-100kb",">100kb")))))[,c(2,3)]
names(dist_table_cut)<-c("dist_label","Freq")
#查看overlab比例：
dist_table_cut[dist_table_cut$dist_label=="has Overlap","Freq"]/sum(dist_table_cut$Freq)


#做图:
png_path<-paste(now_path,"/result.mean_signal.png",sep="")
png(png_path)

p <- ggplot(dist_table_cut,aes(x=dist_label,y=Freq))
p + geom_bar(stat="identity")

p1<-ggplot(dist,aes(dist))
binwidth<-diff(range(dist$dist))/100
p1+geom_density()+xlim(0,100000)

dev.off()

################
#比较心脏的:
rm(list=ls())
dist<-read.table("/home/ding/all_cmd/script/enh_statistics/vista_heart_dist_to_plot",header = F,stringsAsFactors = F)
names(dist)<-c("enh_type","dist")
dist$dist<-abs(dist$dist)
# ggplot(dist,aes(dist))+geom_density()+xlim(0,10^4)

#对距离进行分割:
dist_table<-data.frame(t(table(dist$dist)))[,c(2,3)]
names(dist_table)<-c("dist","Freq")

dist_table_cut<-data.frame(t(table(cut(dist$dist,breaks=c(-1,1,10^3,10^4,10^5,10^6),labels = c("has Overlap","1-1kb","1kb-10kb","10kb-100kb",">100kb")))))[,c(2,3)]
names(dist_table_cut)<-c("dist_label","Freq")
dist_table_cut[dist_table_cut$dist_label=="has Overlap","Freq"]/sum(dist_table_cut$Freq)
dist_table_cut$ratio<-dist_table_cut$Freq/sum(dist_table_cut$Freq)

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/vista_compare_heart.png"
# 使用Cairo生成高清图片：par('din')获得当前
CairoPNG(png_path, width = 6, height = 6, units='in', dpi=600)

p <- ggplot(dist_table_cut,aes(x=dist_label,y=ratio))
p + geom_bar(stat="identity",width=0.6)+
  xlab("distance to vista enhancer")+
  ylab("ratio")+
  annotate("text",x=">100kb",y=0.4,label="Heart",vjust=0.5,size=8,colour="#000080")+
  scale_y_continuous(breaks=seq(0,0.5,0.1),labels=paste(seq(0,50,10),"%",sep=""))+
  theme_bw()+
  theme(axis.text=element_text(size = rel(1.1)),
        axis.text.x=element_text(angle=20,vjust=0.7,hjust=0.5),
        axis.title=element_text(size=rel(1.1)) 
  )
dev.off()

################
#比较脑的:前脑和中脑一起：
rm(list=ls())
dist<-read.table("/home/ding/all_cmd/script/enh_statistics/vista_brain_dist_to_plot",header = F,stringsAsFactors = F)
names(dist)<-c("enh_type","dist")
dist$dist<-abs(dist$dist)
# ggplot(dist,aes(dist))+geom_density()+xlim(0,10^4)

#对距离进行分割:
dist_table<-data.frame(t(table(dist$dist)))[,c(2,3)]
names(dist_table)<-c("dist","Freq")

dist_table_cut<-data.frame(t(table(cut(dist$dist,breaks=c(-1,1,10^3,10^4,10^5,10^6),labels = c("has Overlap","1-1kb","1kb-10kb","10kb-100kb",">100kb")))))[,c(2,3)]
names(dist_table_cut)<-c("dist_label","Freq")
dist_table_cut[dist_table_cut$dist_label=="has Overlap","Freq"]/sum(dist_table_cut$Freq)
dist_table_cut$ratio<-dist_table_cut$Freq/sum(dist_table_cut$Freq)

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/vista_compare_brain.png"
# 使用Cairo生成高清图片：par('din')获得当前
CairoPNG(png_path, width = 6, height = 6, units='in', dpi=600)

p <- ggplot(dist_table_cut,aes(x=dist_label,y=ratio))
p + geom_bar(stat="identity",width=0.6)+
  xlab("distance to vista enhancer")+
  ylab("ratio")+
  annotate("text",x="has Overlap",y=0.6,label="brain",vjust=-0.5,size=8,colour="#000080")+
  scale_y_continuous(breaks=seq(0,0.7,0.1),labels=paste(seq(0,70,10),"%",sep=""))+
  theme_bw()+
  theme(axis.text=element_text(size = rel(1.1)),
        axis.text.x=element_text(angle=20,vjust=0.7,hjust=0.5),
        axis.title=element_text(size=rel(1.1)) 
  )
dev.off()

#########################################
#所有的:
rm(list=ls())
dist<-read.table("/home/ding/all_cmd/script/enh_statistics/dist_to_plot",header = F,stringsAsFactors = F)
names(dist)<-c("enh_type","dist")
dist$dist<-abs(dist$dist)
# ggplot(dist,aes(dist))+geom_density()+xlim(0,10^4)

#对距离进行分割:
dist_table<-data.frame(t(table(dist$dist)))[,c(2,3)]
names(dist_table)<-c("dist","Freq")

dist_table_cut<-data.frame(t(table(cut(dist$dist,breaks=c(-1,1,10^3,10^4,10^5,10^6),labels = c("has Overlap","1-1kb","1kb-10kb","10kb-100kb",">100kb")))))[,c(2,3)]
names(dist_table_cut)<-c("dist_label","Freq")
dist_table_cut[dist_table_cut$dist_label=="has Overlap","Freq"]/sum(dist_table_cut$Freq)
dist_table_cut$ratio<-dist_table_cut$Freq/sum(dist_table_cut$Freq)

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/vista_compare_all.png"
# 使用Cairo生成高清图片：par('din')获得当前
CairoPNG(png_path, width = 6, height = 6, units='in', dpi=600)

p <- ggplot(dist_table_cut,aes(x=dist_label,y=ratio))
p + geom_bar(stat="identity",width=0.6)+
  xlab("distance to vista enhancer")+
  ylab("ratio")+
  # annotate("text",x=">100kb",y=0.4,label="Heart",vjust=0.5,size=8,colour="#000080")+
  scale_y_continuous(breaks=seq(0,0.5,0.1),labels=paste(seq(0,50,10),"%",sep=""))+
  theme_bw()+
  theme(axis.text=element_text(size = rel(1.1)),
        axis.text.x=element_text(angle=20,vjust=0.7,hjust=0.5),
        axis.title=element_text(size=rel(1.1)) 
  )
dev.off()








