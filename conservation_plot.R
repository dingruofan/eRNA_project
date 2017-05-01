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
mean_signal_500<-as.data.frame(read.table("./enh_statistics/conservation_result.mean_signal",header=F,sep="\t",stringsAsFactors = F))
apply(mean_signal_500,2,as.numeric)
names(mean_signal_500)<-c("loc","random","2D-eRNAs Enh","1D-eRNAs Enh","gene TSS","no-eRNAs Enh")


# 平滑：
i<--2000
window_size=10
mean_signal_500_smooth<-data.frame(loc=seq(-2000,2000,1))
while(i<2000){
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"random"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"random"])
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"2D-eRNAs Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"2D-eRNAs Enh"])
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"1D-eRNAs Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"1D-eRNAs Enh"])
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"gene TSS"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"gene TSS"])
  mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"no-eRNAs Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"no-eRNAs Enh"])
  
  i=i+window_size
}
#只能用一个
mean_signal_500_melt<-melt(mean_signal_500_smooth,id=c("loc"))
#只能用一个
# mean_signal_500_melt<-melt(mean_signal_500,id=c("loc"))

names(mean_signal_500_melt)[2]<-c("type")

#排序:
mean_signal_500_melt$type<-factor(mean_signal_500_melt$type,levels = c("2D-eRNAs Enh","1D-eRNAs Enh","no-eRNAs Enh","gene TSS","random"))

#画图:
png_path="/media/ding/000B49000006264C/eRNA_project/figure/conservation_isbid.png"
CairoPNG(png_path, width = 6.2, height = 6, units='in', dpi=600)

ggplot(data=mean_signal_500_melt,aes(x=loc,y=value,colour=type))+
  # geom_smooth()+
  geom_line(size=0.6)+
  # geom_point()+
  xlab("Distance to enhancer center/ TSS of gene (bp)")+
  ylab("PhastCons Score")+
  xlim(c(-1000,1000))+
  scale_colour_manual(limits=c("2D-eRNAs Enh","1D-eRNAs Enh","no-eRNAs Enh","gene TSS","random"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"]),"Gene TSS","random"),values = c("red","gold","green","blue","black"))+
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

# 
# ######################################
# #test:
# 
# rm(list=ls())
# setwd("/media/ding/000B49000006264C/eRNA_project")
# 
# 
# library(ggplot2)
# library(reshape)
# library(Cairo)
# # setwd("/home/ding/all_cmd/script")
# mean_signal_500<-as.data.frame(read.table("./gencode/test.dist",header=F,sep="\t",stringsAsFactors = F))
# apply(mean_signal_500,2,as.numeric)
# names(mean_signal_500)<-c("loc","forward","reverse no stranded","reverse stranded")
# 
# 
# # 平滑：
# # i<--2000
# # window_size=20
# # mean_signal_500_smooth<-data.frame(loc=seq(-2000,2000,1))
# # while(i<2000){
# #   mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"2D-eRNAs Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"2D-eRNAs Enh"])
# #   mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"1D-eRNAs Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"1D-eRNAs Enh"])
# #   mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"gene TSS"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"gene TSS"])
# #   mean_signal_500_smooth[mean_signal_500_smooth$loc %in% seq(i,i+window_size,1),"no-eRNAs Enh"]<-mean(mean_signal_500[mean_signal_500$loc %in% seq(i,i+window_size,1),"no-eRNAs Enh"])
# #   
# #   i=i+window_size
# # }
# # #只能用一个
# # mean_signal_500_melt<-melt(mean_signal_500_smooth,id=c("loc"))
# #只能用一个
# mean_signal_500_melt<-melt(mean_signal_500,id=c("loc"))
# 
# names(mean_signal_500_melt)[2]<-c("type")
# 
# png_path="/media/ding/000B49000006264C/eRNA_project/figure/conservation_gene_TSS.png"
# CairoPNG(png_path, width = par('din')[1], height = par('din')[2], units='in', dpi=600)
# 
# ggplot(data=mean_signal_500_melt,aes(x=loc,y=value,colour=type))+
#   # geom_smooth()+
#   geom_line(size=1,alpha=0.8)+
#   # geom_point()+
#   scale_colour_manual(values=c("blue","green","red"))+
#   xlab("distance to gene TSS(bp)")+
#   ylab("PhastCons Score")+
#   # xlim(c(-1000,1000))+
#   # ggtitle("has strand flag")+
#   # scale_colour_discrete(limits=c("bid_enh","unbid_enh","protein_ss","no_eRNA_enh"),labels=c("2D-eRNAs Enh","1D-eRNAs Enh","gene TSS","no-eRNAs Enh"))+
#   theme(axis.text=element_text(size=rel(1.1)),
#         legend.text=element_text(size=rel(1.)),
#         axis.text=element_text(size = rel(1.1)),
#         axis.title=element_text(size=rel(1.1),),
#         legend.text=element_text(size=rel(1.3)),
#         legend.title=element_blank(),
#         legend.background=element_blank(),
#         legend.key = element_blank(),
#         #legend.position=c(1,1),legend.justification=c(1,1),
#         legend.position="top"
#   )
# dev.off()


#每个组织增强子的conservation boxplot:

rm(list=ls())
setwd("/home/ding/all_cmd/script")

library(ggplot2)
library(reshape)
library(Cairo)
options(bitmapType = "cairo")

phastCons_all<-as.data.frame(read.table("./enh_statistics/conservation/conservation_all.tsv",header=T,sep="\t",stringsAsFactors = F))
phastCons_all$phastCons<-as.numeric(phastCons_all$phastCons)

#标注组织和random不同的颜色:
phastCons_all[phastCons_all$tissue=="random","colours"]<-"grey"
phastCons_all[phastCons_all$tissue!="random","colours"]<-"blue"
# cols<-c("1"="red","0"="green")

#排序：

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/conservation_boxplot.png"
CairoPNG(png_path, width = 5.5, height =5.16 , units='in', dpi=600)

ggplot(phastCons_all,aes(x=reorder(tissue,phastCons,median),y=phastCons))+
  geom_boxplot(aes(fill=colours),show.legend = T,outlier.size = 0,outlier.stroke = 0)+
  ylim(0,0.8)+
  xlab("tissue")+
  ylab("phastCons score")+
  scale_fill_manual(values = c("#FF5151","grey"))+
  # scale_y_continuous(breaks=seq(0,80,20))+
  theme_bw()+
  theme(axis.text=element_text(size=rel(1.1)),
        axis.text.x=element_text(angle = 30,vjust=1,hjust=1),
        legend.text=element_text(size=rel(1.1)),
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

#############
#计算2D-eRNA是否比1D-eRNA要显著保守：
rm(list=ls())
setwd("/home/ding/all_cmd/script/enh_statistics/conservation/conservation_percent")
bid_c<-read.table("bid_enh_conservation",header = T,stringsAsFactors = F)
unbid_c<-read.table("unbid_enh_conservation",header = T,stringsAsFactors = F)
no_eRNA_c<-read.table("no_eRNA_enh_conservation",header = T,stringsAsFactors = F)

all_c<-rbind(bid_c,unbid_c,no_eRNA_c)
t.test(all_c[all_c$type=="bid","phastCons"],all_c[all_c$type=="unbid","phastCons"],alternative = "greater")
t.test(all_c[all_c$type=="unbid","phastCons"],all_c[all_c$type=="no_eRNA","phastCons"],alternative = "greater")



#计算TS是否比UE要显著保守：
rm(list=ls())
setwd("/home/ding/all_cmd/script/enh_statistics/conservation/conservation_percent")
TS_c<-read.table("TS_enh_conservation",header = T,stringsAsFactors = F)
other_c<-read.table("other_enh_conservation",header = T,stringsAsFactors = F)
UE_c<-read.table("UE_enh_conservation",header = T,stringsAsFactors = F)

all_c<-rbind(TS_c,other_c,UE_c)
t.test(all_c[all_c$type=="uni","phastCons"],all_c[all_c$type=="other","phastCons"],alternative = "greater")
t.test(all_c[all_c$type=="other","phastCons"],all_c[all_c$type=="TS","phastCons"],alternative = "greater")
