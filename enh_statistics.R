rm(list=ls())
library(ggplot2)
library(reshape2)
library(Cairo)
library(stringr)
#首先根据聚类信息统计bid_statics和unbid_statistics信息：
enh_count<-read.table("/home/ding/all_cmd/script/enh_statistics/enh_bid_count",header = T,stringsAsFactors = F)
enh_count<-enh_count[,5:14]
tissues<-str_split("Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney"," ")[[1]]
bid_enh<-data.frame(melt(enh_count))
bid_enh<-bid_enh[bid_enh$value>0,]
bid_enh_count<-data.frame(table(bid_enh$variable))
names(bid_enh_count)<-c("tissue","bid_enh_count")

#unbid:
enh_count<-read.table("/home/ding/all_cmd/script/enh_statistics/enh_unbid_count",header = T,stringsAsFactors = F)
enh_count<-enh_count[,5:14]
tissues<-str_split("Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney"," ")[[1]]
unbid_enh<-data.frame(melt(enh_count))
unbid_enh<-unbid_enh[unbid_enh$value>0,]
unbid_enh_count<-data.frame(table(unbid_enh$variable))
names(unbid_enh_count)<-c("tissue","unbid_enh_count")

#all:
enh_count<-read.table("/home/ding/all_cmd/script/enh_statistics/enh_all_count",header = T,stringsAsFactors = F)
enh_count<-enh_count[,5:14]
tissues<-str_split("Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney"," ")[[1]]
all_enh<-data.frame(melt(enh_count))
all_enh<-all_enh[all_enh$value>0,]
all_enh_count<-data.frame(table(all_enh$variable))
names(all_enh_count)<-c("tissue","all_enh_count")

#合并：
bid_statistics<-merge(all_enh_count,bid_enh_count)
bid_statistics$TFBS=NA

unbid_statistics<-merge(all_enh_count,unbid_enh_count)
unbid_statistics$TFBS=NA

#写入文件:
write.table(bid_statistics,"./enh_statistics/bid_statistics",col.names = T,row.names = F,append = F,quote = F)
write.table(unbid_statistics,"./enh_statistics/unbid_statistics",col.names = T,row.names = F,append = F,quote = F)



#############################
rm(list=ls())
setwd("/home/ding/all_cmd/script")



#先读入单向/双向转录enh的数据
bid_statistics<-read.table("./enh_statistics/bid_statistics",stringsAsFactors = F,header = T)
unbid_statistics<-read.table("./enh_statistics/unbid_statistics",stringsAsFactors = F,header = T)

#重新命名列名
names(bid_statistics)<-c("tissue","all_enhancer","enhancer_num","TFBS_num")
names(unbid_statistics)<-c("tissue","all_enhancer","enhancer_num","TFBS_num")

#加一列标注是双向/单向/无eRNA
bid_statistics<-transform(bid_statistics,type="bid_enhancer")
unbid_statistics<-transform(unbid_statistics,type="unbid_enhancer")
noeRNA_statistics<-data.frame(tissue=bid_statistics[,"tissue"],all_enhancer=bid_statistics[,"all_enhancer"],enhancer_num=(bid_statistics[,"all_enhancer"]-bid_statistics[,"enhancer_num"]-unbid_statistics[,"enhancer_num"]),TFBS_num=NA,type="no_eRNA_enhancer")

#合并所有的增强子到一个dataframe
bid_unbid_enh_cbind<-merge(bid_statistics,unbid_statistics,by=c("tissue","all_enhancer"))
bid_unbid_enh_cbind$eRNA_sum<-bid_unbid_enh_cbind$enhancer_num.x+bid_unbid_enh_cbind$enhancer_num.y
sort(bid_unbid_enh_cbind$eRNA_sum)
bid_unbid_enh_cbind$D2_pro<-bid_unbid_enh_cbind$enhancer_num.x/bid_unbid_enh_cbind$all_enhancer
bid_unbid_enh_cbind$eRNA_pro<-bid_unbid_enh_cbind$eRNA_sum/bid_unbid_enh_cbind$all_enhancer

bid_unbid_enh<-rbind(bid_statistics,unbid_statistics,noeRNA_statistics)

#计算每个组织的有eRNA的增强子所占的比例，用来排序和做注释
eRNA_enh_num_ratio<-(bid_unbid_enh[bid_unbid_enh$type=="bid_enhancer","enhancer_num"]+bid_unbid_enh[bid_unbid_enh$type=="unbid_enhancer","enhancer_num"])/bid_unbid_enh[bid_unbid_enh$type=="unbid_enhancer","all_enhancer"]
tissue_order<-bid_unbid_enh$tissue[order(eRNA_enh_num_ratio)]
#将上面排序过的tissue转为因子
bid_unbid_enh$tissue<-factor(bid_unbid_enh$tissue,levels=tissue_order)

#通过annotate添加eRNA和增强子的数目,
num_eRNA<-2*bid_statistics$enhancer_num+unbid_statistics$enhancer_num
num_annot<-data.frame(tissue=append(bid_statistics$tissue,unbid_statistics$tissue),enhancer_num=append(num_eRNA,bid_statistics$all_enhancer))

bid_unbid_enh$type<-factor(bid_unbid_enh$type,levels=c("no_eRNA_enhancer","unbid_enhancer","bid_enhancer"))

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_proportion.png"
CairoPNG(png_path, width = 7.9, height = 6.7, units='in', dpi=600)


ggplot(bid_unbid_enh,aes(x=tissue,y=enhancer_num,fill=type))+
  geom_bar(stat = "identity",position = "fill")+  #stat表示使用原数据,fill表示作为百分比图
  theme_light()+
  coord_flip()+  #反转坐标
  scale_y_continuous(breaks=seq(0,1,0.2),labels=paste(seq(0,100,20),"%",sep=""))+  #设置y轴坐标
  scale_fill_discrete(limits=c("bid_enhancer","unbid_enhancer","no_eRNA_enhancer"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"])))+
  geom_hline(yintercept = 0.5,size=1,color="grey")+
  # xlab("Human tissue")+
  ylab("Percentage of Enhancers")+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1.5)),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=rel(1.5)),
        axis.title.y=element_blank(),
        plot.title=element_text(size=rel(1.5)),
        legend.text=element_text(size=rel(1.5)),
        legend.position="top",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        axis.line=element_blank(),
        axis.ticks.y=element_blank(),
        legend.box.spacing = unit(0.1,"pt"),
        legend.margin = margin(0,0,0.1,0,"pt"),
        # panel.background = element_rect(fill = "red", colour = "grey50"),
        # axis.line.y = element_line( colour = "red"),
        plot.margin=margin(2,53,2,20,unit="pt"),
        panel.border = element_blank()
        # plot.background=element_rect(fill = "black", colour = "grey30")
        )+
  annotate("text",x=1,y=1.001,label=2727,colour="black",size=rel(2.4),hjust=0)+
  annotate("text",x=2,y=1.001,label=12076,colour="black",size=rel(2.4),hjust=0)+
  annotate("text",x=3,y=1.001,label=14193,colour="black",size=rel(2.4),hjust=0)+
  annotate("text",x=4,y=1.001,label=6485,colour="black",size=rel(2.4),hjust=0)+
  annotate("text",x=5,y=1.001,label=11256,colour="black",size=rel(2.4),hjust=0)+
  annotate("text",x=6,y=1.001,label=6476,colour="black",size=rel(2.4),hjust=0)+
  annotate("text",x=7,y=1.001,label=5210,colour="black",size=rel(2.4),hjust=0)+
  annotate("text",x=8,y=1.001,label=11137,colour="black",size=rel(2.4),hjust=0)+
  annotate("text",x=9,y=1.001,label=1374,colour="black",size=rel(2.4),hjust=0)+
  annotate("text",x=10,y=1.001,label=7231,colour="black",size=rel(2.4),hjust=0)


  # annotate("text",x=bid_unbid_enh$tissue,y=bid_unbid_enh$enhancer_num/bid_unbid_enh$all_enhancer,label=bid_unbid_enh$enhancer_num)+

  # ggtitle("Enhacner Num")+

  # annotate(geom="text",x=num_annot$tissue,y=num_annot$enhancer_num,label=num_annot$enhancer_num)
dev.off()

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_number.png"
CairoPNG(png_path, width = 6, height = 6, units='in', dpi=600)

ggplot(bid_unbid_enh,aes(x=reorder(tissue,all_enhancer),y=enhancer_num,fill=type))+
  geom_bar(stat = "identity")+  #stat表示使用原数据,fill表示作为百分比图
  theme_light()+
  coord_flip()+  #反转坐标
  scale_fill_discrete(limits=c("bid_enhancer","unbid_enhancer","no_eRNA_enhancer"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"])))+
  # annotate("text",x=bid_unbid_enh$tissue,y=bid_unbid_enh$enhancer_num/bid_unbid_enh$all_enhancer,label=bid_unbid_enh$enhancer_num)+
  # geom_hline(yintercept = 0.5,size=1,color="grey")+
  xlab("Human tissue")+
  ylab("The number of enhancers")+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1.3),angle=45),
        axis.text.x=element_text(size=rel(1.3)),
        axis.title.x=element_text(size=rel(1.3)),
        axis.title.y=element_blank(),
        plot.title=element_text(size=rel(1)),
        legend.position="top",
        legend.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.box.spacing = unit(1,"mm")
  )
# annotate(geom="text",x=num_annot$tissue,y=num_annot$enhancer_num,label=num_annot$enhancer_num)
dev.off()

################################################################
#做图：增强子在各个组织中的分布数目图
#包括三条线：bid enh/unbid enh/noeRNA enh/all enh增强子的直方图/
rm(list=ls())
setwd("/home/ding/all_cmd/script")

library(ggplot2)
library(Cairo)
enh_bid_count<-read.table("./enh_statistics/enh_bid_count",header = T,stringsAsFactors = F)
enh_unbid_count<-read.table("./enh_statistics/enh_unbid_count",header = T,stringsAsFactors = F)
enh_no_eRNA_count<-read.table("./enh_statistics/enh_no_eRNA_count",header = T,stringsAsFactors = F)

enh_bid_count$type="bid"
enh_unbid_count$type="unbid"
enh_no_eRNA_count$type="no_eRNA"

enh_all_count<-rbind(enh_bid_count,enh_unbid_count,enh_no_eRNA_count)
write.table(enh_all_count,"./enh_statistics/enh_all_count",append = F,quote = F,sep="\t",row.names = F,col.names = T)

enh_bid_table<-data.frame(t(table(enh_bid_count$tissue_statistics)))
enh_bid_table$Var1<-"bid"
enh_unbid_table<-data.frame(t(table(enh_unbid_count$tissue_statistics)))
enh_unbid_table$Var1<-"unbid"
enh_unbid_table<-rbind(enh_unbid_table,data.frame(Var1="unbid",Var2="9",Freq=0))
enh_unbid_table<-rbind(enh_unbid_table,data.frame(Var1="unbid",Var2="10",Freq=0))
enh_no_eRNA_table<-data.frame(t(table(enh_no_eRNA_count$tissue_statistics)))
enh_no_eRNA_table$Var1<-"no_eRNA"

enh_all_table<-rbind(enh_bid_table,enh_unbid_table,enh_no_eRNA_table)
enh_all_table$Var2=factor(enh_all_table$Var2,levels=1:10)

enh_all_table[enh_all_table$Var1=="bid","Var1"]<-"2D-eRNAs Enh"
enh_all_table[enh_all_table$Var1=="unbid","Var1"]<-"1D-eRNAs Enh"
enh_all_table[enh_all_table$Var1=="no_eRNA","Var1"]<-"no-eRNAs Enh"

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_tissue_num_propotion.png"
CairoPNG(png_path, width = 6, height = 6, units='in', dpi=600)

enh_all_table$Var1<-factor(enh_all_table$Var1,levels=c("no-eRNAs Enh","1D-eRNAs Enh","2D-eRNAs Enh"))


ggplot(enh_all_table,aes(x=Var2,y=Freq,fill=Var1))+
  geom_bar(stat = "identity",position  = "fill")+ 
  scale_fill_discrete(limits=c("2D-eRNAs Enh","1D-eRNAs Enh","no-eRNAs Enh"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"])))+
  scale_y_continuous(breaks=seq(0,1,0.2),labels = paste(seq(0,100,20),"%",sep=""))+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1.3)),
        axis.text.x=element_text(size=rel(1.3)),
        axis.title.x=element_text(size=rel(1.3)),
        axis.title.y=element_text(size=rel(1.3)),
        legend.position="top",
        legend.title=element_blank(),
        legend.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.spacing = unit(1,"mm"),
        legend.background=element_blank(),
        legend.key = element_blank())+
  xlab("The number of tissues")+
  ylab("Percentage of enhancers")
  # ggtitle("Enhacner distribute in tissues")+theme(plot.title=element_text(size=rel(2)))
dev.off()


png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_tissue_num.png"
CairoPNG(png_path, width = 6, height = 6, units='in', dpi=600)

enh_all_table_1<-enh_all_table
enh_all_table_1[enh_all_table_1$Freq==0,"Freq"]<-1
enh_all_table_1[enh_all_table_1$Freq==1,"Freq"]<-2
enh_all_table_1$Var1<-factor(enh_all_table_1$Var1,levels=c("2D-eRNAs Enh","1D-eRNAs Enh","no-eRNAs Enh"))

ggplot(enh_all_table_1,aes(x=Var2,y=Freq,fill=Var1))+
  # geom_bar(stat = "identity",position  = "fill")+ 
  geom_bar(stat = "identity",position="dodge",width=0.8)+ 
  scale_y_log10()+
  scale_fill_discrete(limits=c("2D-eRNAs Enh","1D-eRNAs Enh","no-eRNAs Enh"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"])))+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1.3)),
        axis.text.x=element_text(size=rel(1.3)),
        axis.title.x=element_text(size=rel(1.3)),
        axis.title.y=element_text(size=rel(1.3)),
        legend.text = element_text(size=rel(1.3)),
        legend.position="top",
        legend.title=element_blank(),
        legend.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.spacing = unit(1,"mm"),
        legend.background=element_blank(),
        legend.key = element_blank())+
  xlab("The number of tissues ")+
  ylab("The number of Enhancers ")

# ggtitle("Enhacner distribute in tissues")+theme(plot.title=element_text(size=rel(2)))
dev.off()





png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_tissue_num_isbid.png"
CairoPNG(png_path, width = 6, height = 6, units='in', dpi=600)

b<-enh_all_table
#转换为因子并重新定义label:
b$Var1<-factor(b$Var1,levels=c("2D-eRNAs Enh","1D-eRNAs Enh","no-eRNAs Enh"),labels=c("Enh['2D-eRNA']","Enh['1D-eRNA']","Enh['no-eRNA']"))
#将Freq数值为1的替换为2,避免log10为0:为0的删除。
b[b$Freq==1,"Freq"]=2
b<-b[b$Freq>0,]

ggplot(b,aes(x=Var2,y=Freq,fill=Var1))+ 
  geom_bar(stat = "identity")+
  scale_y_log10()+
  facet_wrap(~Var1,scales="free",labeller = label_parsed)+
  guides(fill=FALSE)+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1.0)),
        axis.text.x=element_text(size=rel(1.0)),
        axis.title.x=element_text(size=rel(1.1)),
        axis.title.y=element_text(size=rel(1.1)),
        legend.position="top",
        legend.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        strip.text.x=element_text(size=rel(1.1))
        )+
  xlab("Human tissue Number")+
  ylab("Enhancer number")
  # ggtitle("Enhacner Num")+theme(plot.title=element_text(size=rel(2)),legend.position="top",legend.title=element_blank())
dev.off()

