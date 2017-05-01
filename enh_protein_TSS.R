rm(list=ls())
library(ggplot2)
library(Cairo)
setwd("/home/ding/all_cmd/script/enh_statistics/")

enh_dist<-read.table("./enh_protein_TSS_dist",header = F,stringsAsFactors = F)
names(enh_dist)<-c("type","dist")
# enh_dist<-subset(enh_dist,abs(dist)<=300000)

binsize<-diff(range(enh_dist$dist))/1000
# ggplot(enh_dist,aes(x=dist))+geom_histogram(binwidth=binsize,aes(fill=type),alpha=0.5)+xlim(-500000,500000)
table_bid<-data.frame(t(table(cut(enh_dist[enh_dist$type=="bid","dist"],breaks=c(c(-6,-1)*10^6,-5*10^5,-10^5,-10^4,-10^3,0,10^3,10^4,10^5,5*10^5,c(1,6)*10^6),labels=c("<-1000kb","-1000kb~-500kb","-500kb~-100kb","-100kb~-10kb","-10kb~-1kb","-1kb~-0","0~1kb","1kb~10kb","10kb~100kb","100kb~500kb","500kb~1000kb",">1000kb")))),type="bid")[c(2:4)]
table_unbid<-data.frame(t(table(cut(enh_dist[enh_dist$type=="unbid","dist"],breaks=c(c(-6,-1)*10^6,-5*10^5,-10^5,-10^4,-10^3,0,10^3,10^4,10^5,5*10^5,c(1,6)*10^6),labels=c("<-1000kb","-1000kb~-500kb","-500kb~-100kb","-100kb~-10kb","-10kb~-1kb","-1kb~-0","0~1kb","1kb~10kb","10kb~100kb","100kb~500kb","500kb~1000kb",">1000kb")))),type="unbid")[c(2:4)]
table_no_eRNA<-data.frame(t(table(cut(enh_dist[enh_dist$type=="no_eRNA","dist"],breaks=c(c(-6,-1)*10^6,-5*10^5,-10^5,-10^4,-10^3,0,10^3,10^4,10^5,5*10^5,c(1,6)*10^6),labels=c("<-1000kb","-1000kb~-500kb","-500kb~-100kb","-100kb~-10kb","-10kb~-1kb","-1kb~-0","0~1kb","1kb~10kb","10kb~100kb","100kb~500kb","500kb~1000kb",">1000kb")))),type="no_eRNA")[c(2:4)]

table_all<-rbind(table_bid,table_unbid,table_no_eRNA)
names(table_all)<-c("dist","Freq","type")

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_TSS_dist_density.png"
CairoPNG(png_path, width = 6, height = 5.5, units='in', dpi=600)

# ggplot(enh_dist,aes(x=type,y=dist))+geom_line(stat="identity")+expand_limits(y=0)
ggplot(enh_dist,aes(x=dist))+
  geom_density(aes(color=type),alpha=0.8)+
  scale_x_continuous(limits = c(-500000,500000),breaks=seq(-500000,500000,100000),labels = paste(seq(-500,500,100),sep=""))+
  xlab("distance to proximate gene (kb)")+
  ylab("density")+
  scale_colour_discrete(limits=c("bid","unbid","no_eRNA"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"])))+
  theme_bw()+
  theme(
        axis.text=element_text(size=rel(0.9)),
        axis.title=element_text(size=rel(1.0)),
        legend.text=element_text(size=rel(1.0)),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.margin=unit(0,"mm"),
        #legend.position=c(1,1),legend.justification=c(1,1),
        legend.position="top"
  )

dev.off()


png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_TSS_dist_density_all_100kb.png"
CairoPNG(png_path, width = 7.2, height = 6.6, units='in', dpi=600)

# ggplot(enh_dist,aes(x=type,y=dist))+geom_line(stat="identity")+expand_limits(y=0)
ggplot(enh_dist,aes(x=dist))+
  # geom_histogram(binwidth =diff(range(enh_dist$dist))/4022 ,alpha=0.3)+
  geom_density(alpha=0.8,adjust=1/1000000)+
  scale_x_continuous(limits = c(-100000,100000),breaks=c(seq(-100000,100000,20000)),labels = paste(c(seq(-100,100,20)),sep=""))+
  # geom_vline(xintercept=-1000)+
  xlab("The Distance between enhancers and their adjacent genes(kb)")+
  ylab("Density")+
  scale_colour_discrete(limits=c("bid","unbid","no_eRNA"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"])))+
  theme_bw()+
  theme(
    axis.text=element_text(size=rel(1.2)),
    axis.title=element_text(size=rel(1.3)),
    legend.text=element_text(size=rel(1.3)),
    legend.title=element_blank(),
    legend.background=element_blank(),
    legend.key = element_blank(),
    #legend.position=c(1,1),legend.justification=c(1,1),
    legend.position="top"
  )

dev.off()


########做点图:
enh_dist_table<-as.data.frame(table(enh_dist$dist),stringsAsFactors = F)
names(enh_dist_table)<-c("loc","Freq")
enh_dist_table$loc<-as.numeric(enh_dist_table$loc)
enh_dist_table$Freq<-as.numeric(enh_dist_table$Freq)
#选择上下游500kb的做图:
enh_dist_table<-subset(enh_dist_table,loc<(500000)&loc>(-500000))

# 补全loc:
enh_dista_table_com<-data.frame(loc=seq(range(enh_dist_table$loc)[1],range(enh_dist_table$loc)[2],1))
enh_dista_table_com_all<-merge(enh_dista_table_com,enh_dist_table,by="loc",all.x = T)
enh_dista_table_com_all[is.na(enh_dista_table_com_all$Freq),"Freq"]<-0

#使用1kb进行平滑:

# 平滑：
i<--500000
window_size=1000
mean_signal_500_smooth<-data.frame(loc=seq(-500000,500000,1))
while(i<500000){
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"Freq"]<-mean(enh_dista_table_com_all[enh_dista_table_com_all$loc %in% seq(i,i+window_size,1),"Freq"])
  i=i+window_size
  print(i)
}


ggplot(enh_dista_table_com_all,aes(x=loc,y=Freq))+
  geom_bar(stat="identity",alpha=0.8)+
  scale_x_continuous(limits = c(-50000,50000))+
  xlab("distance to proximate gene (kb)")+
  ylab("density")+
  scale_colour_discrete(limits=c("bid","unbid","no_eRNA"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"])))+
  theme_bw()+
  theme(
    axis.text=element_text(size=rel(0.9)),
    axis.title=element_text(size=rel(1.0)),
    legend.text=element_text(size=rel(1.0)),
    legend.title=element_blank(),
    legend.background=element_blank(),
    legend.key = element_blank(),
    legend.margin=unit(0,"mm"),
    #legend.position=c(1,1),legend.justification=c(1,1),
    legend.position="top"
  )



########################################################

rm(list=ls())
library(ggplot2)
library(Cairo)
setwd("/home/ding/all_cmd/script/enh_statistics/")

enh_dist<-read.table("./enh_protein_TSS_dist2",header = F,stringsAsFactors = F)
names(enh_dist)<-c("type","dist")
# enh_dist<-subset(enh_dist,abs(dist)<=300000)

binsize<-diff(range(enh_dist$dist))/1000
# ggplot(enh_dist,aes(x=dist))+geom_histogram(binwidth=binsize,aes(fill=type),alpha=0.5)+xlim(-500000,500000)
table_bid<-data.frame(t(table(cut(enh_dist[enh_dist$type=="spe","dist"],breaks=c(c(-6,-1)*10^6,-5*10^5,-10^5,-10^4,-10^3,0,10^3,10^4,10^5,5*10^5,c(1,6)*10^6),labels=c("<-1000kb","-1000kb~-500kb","-500kb~-100kb","-100kb~-10kb","-10kb~-1kb","-1kb~-0","0~1kb","1kb~10kb","10kb~100kb","100kb~500kb","500kb~1000kb",">1000kb")))),type="bid")[c(2:4)]
table_unbid<-data.frame(t(table(cut(enh_dist[enh_dist$type=="oth","dist"],breaks=c(c(-6,-1)*10^6,-5*10^5,-10^5,-10^4,-10^3,0,10^3,10^4,10^5,5*10^5,c(1,6)*10^6),labels=c("<-1000kb","-1000kb~-500kb","-500kb~-100kb","-100kb~-10kb","-10kb~-1kb","-1kb~-0","0~1kb","1kb~10kb","10kb~100kb","100kb~500kb","500kb~1000kb",">1000kb")))),type="unbid")[c(2:4)]
table_no_eRNA<-data.frame(t(table(cut(enh_dist[enh_dist$type=="uni","dist"],breaks=c(c(-6,-1)*10^6,-5*10^5,-10^5,-10^4,-10^3,0,10^3,10^4,10^5,5*10^5,c(1,6)*10^6),labels=c("<-1000kb","-1000kb~-500kb","-500kb~-100kb","-100kb~-10kb","-10kb~-1kb","-1kb~-0","0~1kb","1kb~10kb","10kb~100kb","100kb~500kb","500kb~1000kb",">1000kb")))),type="no_eRNA")[c(2:4)]

table_all<-rbind(table_bid,table_unbid,table_no_eRNA)
names(table_all)<-c("dist","Freq","type")

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_TSS_dist_density2.png"
CairoPNG(png_path, width = 6, height = 5.5, units='in', dpi=600)

# ggplot(enh_dist,aes(x=type,y=dist))+geom_line(stat="identity")+expand_limits(y=0)
ggplot(enh_dist,aes(x=dist))+
  geom_density(aes(color=type),alpha=0.8)+
  scale_x_continuous(limits = c(-500000,500000),breaks=seq(-500000,500000,100000),labels = paste(seq(-500,500,100),sep=""))+
  xlab("distance to proximate gene (kb)")+
  ylab("density")+
  scale_colour_discrete(limits=c("spe","oth","uni"),labels=c(expression(Enh["TS"]),expression(Enh["Oth"]),expression(Enh["UE"])))+
  theme_bw()+
  theme(axis.text=element_text(size=rel(0.9)),
        legend.text=element_text(size=rel(1.0)),
        axis.title=element_text(size=rel(1.0)),
        legend.text=element_text(size=rel(1.0)),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.margin=unit(0,"mm"),
        #legend.position=c(1,1),legend.justification=c(1,1),
        legend.position="top"
  )
dev.off()





  