#从shell读入工作目录
rm(list=ls())
library(ggplot2)
library(reshape)
library(Cairo)
setwd("/home/ding/all_cmd/script/enh_statistics")

dist_iseRNA<-read.table("./TF_enh_dist_iseRNA_type",header = F,stringsAsFactors = F)
names(dist_iseRNA)<-c("TF","type","dist")
dist_isspe<-read.table("./TF_enh_dist_isspe_type",header = F,stringsAsFactors = F)
names(dist_isspe)<-c("TF","type","dist")

binwidth_1<-diff(range(dist_iseRNA$dist))
dist_iseRNA$type<-factor(dist_iseRNA$type,levels = c("bid","unbid","no_eRNA"))

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/TF_enh_dist_iseRNA.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 6, height = 6, units='in', dpi=600)

ggplot(dist_iseRNA,aes(x=dist,y=..density..))+
  theme_bw()+
  geom_density(position = "identity",aes(color=type),alpha=1)+
  xlab("distance to enhancer center (bp)")+
  scale_x_continuous(breaks=c(-2000,-400,0,400,2000),labels=c(-2000,-400,0,400,2000),limits = c(-2400,2400))+
  geom_vline(xintercept = -400,color="#AAA000",linetype="dotted")+
  geom_vline(xintercept = 400,color="#AAA000",linetype="dotted")+
  scale_color_discrete(breaks=c("bid","unbid","no_eRNA"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"])))+
  theme_bw()+
  theme(
        axis.text=element_text(size=rel(1.3)),
        # panel.border=element_blank(),
        legend.title=element_blank(),
        legend.margin=margin(0,0,0,0,"mm"),
        legend.position="top",
        legend.text=element_text(size=rel(1.3)),
        legend.box.spacing = unit(1,"mm")
  )
dev.off()


######################################
#不分类的结果:

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/TF_enh_dist_all.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 6, height = 6, units='in', dpi=600)

# dist_iseRNA_1<-dist_iseRNA[dist_iseRNA$dist<=2000&dist_iseRNA$dist>=-2000,]

# binwidth_1<-diff(range(dist_iseRNA$dist))/0

ggplot(dist_iseRNA,aes(x=dist,y=..density..))+
  theme_bw()+
  geom_density(position = "identity",alpha=1)+
  # geom_histogram(position = "identity",binwidth = 50,alpha=0.3)+
  xlab("The distance between TFBS and enhancer center (bp)")+
  ylab("Density")+
  scale_x_continuous(breaks=c(-10000,-5000,-500,0,500,5000,10000),labels=c(-10000,-5000,-500,0,500,5000,10000))+
  geom_vline(xintercept = -500,color="#000000",linetype="dotted")+
  geom_vline(xintercept = 500,color="#000000",linetype="dotted")+
  theme_bw()+
  theme(
    axis.text.x=element_text(size=rel(1.1),angle = 45,hjust = 1,vjust=1),
    axis.title = element_text(size=rel(1.3)),
    axis.text.y = element_text(size=rel(1.3)),
    # panel.border=element_blank(),
    legend.title=element_blank(),
    legend.margin=margin(0,0,0,0,"mm"),
    legend.position="top",
    legend.text=element_text(size=rel(1.3)),
    legend.box.spacing = unit(1,"mm")
  )

dev.off()

############
#直接做点图bid:
library(reshape)
dist_iseRNA_table<-as.data.frame.ts(table(dist_iseRNA$dist,dist_iseRNA$type))
dist_iseRNA_table$location<-rownames(dist_iseRNA_table)

# 平滑：
i<--9000
window_size=10
dist_iseRNA_table_smoth<-data.frame(location=seq(-9000,9000,1))
while(i<9000){
  dist_iseRNA_table_smoth[dist_iseRNA_table_smoth$location %in% seq(i,i+window_size,1),"bid"]<-mean(dist_iseRNA_table[dist_iseRNA_table$location %in% seq(i,i+window_size,1),"bid"])
  dist_iseRNA_table_smoth[dist_iseRNA_table_smoth$location %in% seq(i,i+window_size,1),"unbid"]<-mean(dist_iseRNA_table[dist_iseRNA_table$location %in% seq(i,i+window_size,1),"unbid"])
  dist_iseRNA_table_smoth[dist_iseRNA_table_smoth$location %in% seq(i,i+window_size,1),"no_eRNA"]<-mean(dist_iseRNA_table[dist_iseRNA_table$location %in% seq(i,i+window_size,1),"no_eRNA"])
  
  i=i+window_size
}

dist_iseRNA_table_smoth_melt<-melt(dist_iseRNA_table_smoth,id=c("location"))
names(dist_iseRNA_table_smoth_melt)<-c("location","type","Freq")

ggplot(dist_iseRNA_table_smoth_melt,aes(x=location,y=Freq))+geom_line(aes(colour=type))+ylim(0,20)+theme_bw()+xlim(-2000,2000)




############
#直接做点图spe:
library(reshape)
dist_isspe_table<-as.data.frame.ts(table(dist_isspe$dist,dist_isspe$type))
dist_isspe_table$location<-rownames(dist_isspe_table)

# 平滑：
i<--9000
window_size=10
dist_isspe_table_smoth<-data.frame(location=seq(-9000,9000,1))
while(i<9000){
  dist_isspe_table_smoth[dist_isspe_table_smoth$location %in% seq(i,i+window_size,1),"other"]<-mean(dist_isspe_table[dist_isspe_table$location %in% seq(i,i+window_size,1),"other"])
  dist_isspe_table_smoth[dist_isspe_table_smoth$location %in% seq(i,i+window_size,1),"spe"]<-mean(dist_isspe_table[dist_isspe_table$location %in% seq(i,i+window_size,1),"spe"])
  dist_isspe_table_smoth[dist_isspe_table_smoth$location %in% seq(i,i+window_size,1),"uni"]<-mean(dist_isspe_table[dist_isspe_table$location %in% seq(i,i+window_size,1),"uni"])
  
  i=i+window_size
}

dist_isspe_table_smoth<-melt(dist_iseRNA_table_smoth,id=c("location"))
names(dist_isspe_table_smoth)<-c("location","type","Freq")

ggplot(dist_isspe_table_smoth,aes(x=location,y=Freq))+geom_line(aes(colour=type))+ylim(0,20)+theme_bw()+xlim(-2000,2000)

################################
#其实不用分类:
ggplot(dist_isspe_table_smoth,aes(x=location,y=Freq))+geom_density()+ylim(0,20)+theme_bw()+xlim(-2000,2000)



##########################################
png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/TF_enh_dist_isspe.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 6, height = 6, units='in', dpi=600)

ggplot(dist_isspe,aes(x=dist,y=..density..))+
  theme_bw()+
  geom_density(position = "identity",aes(color=type),alpha=1)+
  xlab("distance to enhancer center (bp)")+
  scale_x_continuous(breaks=c(-2000,-400,0,400,2000),labels=c(-2000,-400,0,400,2000),limits = c(-2400,2400))+
  geom_vline(xintercept = -400,color="#AAA000",linetype="dotted")+
  geom_vline(xintercept = 400,color="#AAA000",linetype="dotted")+
  scale_color_discrete(breaks=c("spe","other","uni"),labels=c(expression(Enh["TS"]),expression(Enh["Oth"]),expression(Enh["UE"])))+
  theme_bw()+
  theme(
    axis.text.x=element_text(size=rel(1.1)),
    # panel.border=element_blank(),
    legend.title=element_blank(),
    legend.margin=unit(0.1,"mm"),
    legend.position="top",
    legend.text=element_text(size=rel(1.1))
  )
dev.off()

