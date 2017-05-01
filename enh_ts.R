#从shell读入工作目录
rm(list=ls())
library(ggplot2)
library(Cairo)
setwd("/home/ding/all_cmd/script")

#必须先运行enh_statistics.R生成enh_all_count文件
histone_matrix<-read.table("./enh_statistics/enh_all_count",header = T,stringsAsFactors = F)
#根据公式计算Pts:
tissues<-names(histone_matrix)[c(-1:-4,-ncol(histone_matrix),-ncol(histone_matrix)+1)]
rowsum_sig<-as.matrix(rowSums(histone_matrix[,tissues,drop=F]),ncol=1)

for(i in tissues){
  histone_matrix[,i]<-histone_matrix[,i,drop=F]/rowsum_sig
}

#计算香农熵Hs作为组织特异性分值:
for(i in tissues){
  histone_matrix[,i]<-log(histone_matrix[,i,drop=F]^histone_matrix[,i,drop=F],2)
}
ts_matrix<-subset(histone_matrix,select = c("chr","chr_start","chr_end","enhancer","type"))
ts_matrix$Hs<-(-1)*rowSums(histone_matrix[,tissues,drop=F])


quantile(ts_matrix[ts_matrix$type=="bid","Hs"],probs = seq(0,1,0.2))
quantile(ts_matrix[ts_matrix$type=="unbid","Hs"],probs = seq(0,1,0.2))
quantile(ts_matrix[ts_matrix$type=="no_eRNA","Hs"],probs = seq(0,1,0.2))

#取80%分位数作为阈值，小于阈值的是特异性enh:
#小于1的是spe，取1-2的是other,大于2的是uni:
# specific_enh_score<-quantile(ts_matrix[,"Hs"],probs = seq(0,1,0.2))[5]
specific_enh_score<-0.5

ts_matrix$type<-factor(ts_matrix$type,levels=c("no_eRNA","unbid","bid"))

#做图:
png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/ts_enh_density.png"
CairoPNG(png_path, width = 6.3, height = 6.3, units='in', dpi=600)

ggplot(ts_matrix,aes(Hs))+
  geom_density(aes(fill=type),alpha=0.5)+
  # geom_vline(xintercept=specific_enh_score,colour="#AAA000", linetype="dashed")+
  # geom_vline(xintercept=1.8,colour="#AAA000", linetype="dashed")+
  scale_x_continuous(breaks=round(as.vector(sort(c(seq(0,3,1),specific_enh_score))),0))+
  # scale_fill_discrete()+
  scale_fill_manual(breaks=c("bid","unbid","no_eRNA"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"])),
                    values=c("#7CAE00","#00BFC4","#F8766D"))+
  #values = c( "#ff0000", "#00ff00","#eec900"
  xlab("Tissue specificity index of enhancer")+
  ylab("Density")+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1.3)),
        axis.text.x=element_text(size=rel(1.3)),
        axis.title.x=element_text(size=rel(1.3)),
        axis.title.y=element_text(size=rel(1.3)),
        plot.title=element_text(size=rel(1.3)),
        legend.text = element_text(size=rel(1.3)),
        legend.position="top",
        legend.margin=margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,0,0,0,"pt"),
        legend.box.spacing = unit(0.2,"cm"),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank()
  )

dev.off()

ts_matrix[which(ts_matrix$Hs<specific_enh_score),"is_spe"]<-"spe"
ts_matrix[which(ts_matrix$Hs<2&ts_matrix$Hs>=specific_enh_score),"is_spe"]<-"other"
ts_matrix[which(ts_matrix$Hs>=2),"is_spe"]<-"uni"

#将spe/uni/other enh写入不同文件：
enh_spe<-ts_matrix[ts_matrix$is_spe=="spe",c(1,2,3,4)]
enh_uni<-ts_matrix[ts_matrix$is_spe=="uni",c(1,2,3,4)]
enh_other<-ts_matrix[ts_matrix$is_spe=="other",c(1,2,3,4)]

write.table(enh_spe,"./enh_statistics/enh_spe_cluster.bed",append=F,quote = F,sep="\t",col.names = F,row.names = F)
write.table(enh_uni,"./enh_statistics/enh_uni_cluster.bed",append=F,quote = F,sep="\t",col.names = F,row.names = F)
write.table(enh_other,"./enh_statistics/enh_other_cluster.bed",append=F,quote = F,sep="\t",col.names = F,row.names = F)

#将所有enh信息写入enh_spe:
write.table(ts_matrix,"./enh_statistics/enh_spe",append=F,quote = F,sep="\t",col.names = T,row.names = F)
ts_matrix$type<-factor(ts_matrix$type,levels = c("bid","unbid","no_eRNA"))
ts_matrix$is_spe<-factor(ts_matrix$is_spe,levels = c("spe","other","uni"))

#做图:
png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/ts_enh_num.png"
CairoPNG(png_path, width = 6, height = 6, units='in', dpi=600)

ggplot(ts_matrix,aes(x=is_spe))+
  geom_bar(aes(fill=type),position = "dodge",width=0.7)+
  # xlab("Tissue specificity index")+
  ylab("enhancer number")+
  scale_fill_discrete(breaks=c("bid","unbid","no_eRNA"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"])))+
  scale_y_log10(limits=c(1,10^5),breaks=10^(0:5),labels=c(1,10^(1:4),as.integer(100000)))+
  scale_x_discrete(breaks=c("spe","other","uni"),labels=c(expression(Enh["TS"]),expression(Enh["Oth"]),expression(Enh["UE"])))+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1.0)),
        axis.text.x=element_text(size=rel(1.1)),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=rel(1.1)),
        plot.title=element_text(size=rel(1)),
        legend.position="top",
        legend.margin=unit(0.1,"mm"),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank()
  )

dev.off()

##############################################
##############################################
#enh注释：添加到/home/ding/all_cmd/script/enh_statistics/enh_all_count_spe
rm(list=ls())
library(reshape2)
library(ggplot2)
library(Cairo)
setwd("/home/ding/all_cmd/script")
#读入数据:
histone_matrix<-read.table("./enh_statistics/enh_all_count",header = T,stringsAsFactors = F)
ts_matrix<-read.table("./enh_statistics/enh_spe",header = T,stringsAsFactors = F)
#去掉Hs列：
ts_matrix<-ts_matrix[c(1,2,3,4,5,7)]
#进行merge:
enh_all_count_spe<-merge(histone_matrix,ts_matrix,by=c("chr","chr_start","chr_end","enhancer","type"))
#调整列的顺序:
enh_all_count_spe<-enh_all_count_spe[c(1:4,6:16,5,17)]
#写入文件:
write.table(enh_all_count_spe,"./enh_statistics/enh_all_count_spe",append=F,quote = F,sep="\t",col.names = T,row.names = F)
#写入文件:不包括组织信息：
enh_all_count_spe_no_tissue<-enh_all_count_spe[,c("chr","chr_start","chr_end","enhancer","type","is_spe")]
write.table(enh_all_count_spe_no_tissue,"./enh_statistics/enh_all",append=F,quote = F,sep="\t",col.names = T,row.names = F)


#宽格式变长格式后写入文件:
enh_all_count_spe_melt<-melt(enh_all_count_spe,id.vars = c("chr","chr_start","chr_end","enhancer","tissue_statistics","type","is_spe"))
names(enh_all_count_spe_melt)[8]<-"tissue"
enh_all_count_spe_melt<-enh_all_count_spe_melt[enh_all_count_spe_melt$value==1,]
enh_all_count_spe_melt$value<-NULL
#写入文件:
write.table(enh_all_count_spe_melt,"./enh_statistics/enh_all_count_spe_melt",append=F,quote = F,sep="\t",col.names = T,row.names = F)

#做图:
enh_all_count_spe_melt_table<-data.frame(table(enh_all_count_spe_melt$tissue,enh_all_count_spe_melt$is_spe))
names(enh_all_count_spe_melt_table)<-c("tissue","is_spe","Freq")
tissue_order_factor_level<-enh_all_count_spe_melt_table[order(enh_all_count_spe_melt_table[enh_all_count_spe_melt_table$is_spe=="spe","Freq"]),"tissue"]
enh_all_count_spe_melt_table$tissue<-factor(enh_all_count_spe_melt_table$tissue,levels = tissue_order_factor_level)

#根据是否特异做增强子的统计图。
#数量图:
png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_number_isspe.png"
CairoPNG(png_path, width = 6, height = 6, units='in', dpi=600)

ggplot(enh_all_count_spe_melt_table,aes(x=reorder(tissue,Freq),y=Freq))+
  geom_bar(stat = "identity",aes(fill=is_spe))+
  coord_flip()+
  ylab("The number of enhancers")+
  
  scale_fill_discrete(c=100,breaks=c("spe","other","uni"),labels=c(expression( Enh["TS"]),expression(Enh["Oth"]),expression( Enh["UE"])))+
  scale_color_manual(values = c("red", "green","blue"))+
  # scale_y_log10(limits=c(1,10^5),breaks=10^(0:5),labels=c(1,10^(1:4),as.integer(100000)))+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1.3),angle = 30),
        axis.text.x=element_text(size=rel(1.3)),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=rel(1.3)),
        plot.title=element_text(size=rel(1.3)),
        legend.position="top",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,0, 0),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.box.spacing = unit(1,"mm")
  )

dev.off()

#根据是否特异做增强子的统计图。
#百分比图:

#做图:


enh_all_count_spe_melt_table<-data.frame(table(enh_all_count_spe_melt$tissue,enh_all_count_spe_melt$is_spe))
names(enh_all_count_spe_melt_table)<-c("tissue","is_spe","Freq")
enh_all_count_spe_melt_table$tissue<-factor(enh_all_count_spe_melt_table$tissue,levels = tissue_order_factor_level)

#根据是否特异做增强子的统计图。
#数量图:
png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_propotion_isspe.png"
CairoPNG(png_path, width = 6.2, height = 6, units='in', dpi=600)

#排序：
enh_all_count_spe_melt_table$is_spe<-factor(enh_all_count_spe_melt_table$is_spe,levels=c("uni","other","spe"))
#prop:
enh_all_count_spe_melt_table_prop<-data.frame(prop.table(table(enh_all_count_spe_melt$tissue,enh_all_count_spe_melt$is_spe),1))
names(enh_all_count_spe_melt_table_prop)<-c("tissue","is_spe","prop")

enh_all_count_spe_melt_table_prop_spe<-subset(enh_all_count_spe_melt_table_prop,is_spe=="spe")
tissue_ordered<-levels(reorder(enh_all_count_spe_melt_table_prop_spe$tissue,enh_all_count_spe_melt_table_prop_spe$prop))

enh_all_count_spe_melt_table$tissue<-factor(enh_all_count_spe_melt_table$tissue,levels = tissue_ordered)

ggplot(enh_all_count_spe_melt_table,aes(x=tissue,y=Freq))+
  geom_bar(stat = "identity",aes(fill=is_spe),position = "fill")+
  coord_flip()+
  ylab("Percentage of Enhancers")+
  scale_fill_discrete(breaks=c("spe","other","uni"),labels=c(expression( Enh["TS"]),expression(Enh["Oth"]),expression( Enh["UE"])))+
  scale_y_continuous(breaks=seq(0,1,0.2),labels=paste(seq(0,100,20),"%",sep=""))+
  # scale_y_log10(limits=c(1,10^5),breaks=10^(0:5),labels=c(1,10^(1:4),as.integer(100000)))+
  theme_classic()+
  theme(axis.text.x=element_text(size=rel(1.3)),
        axis.text.y=element_text(size=rel(1.3),vjust = 1,hjust = 1),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=rel(1.3)),
        legend.text = element_text(size=rel(1.3)),
        legend.position="top",
        legend.margin = margin(0,0,0,0),
        legend.box.spacing =margin(0,0,0, 0),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        plot.margin=margin(2,33,2,2,unit="pt"),
        axis.line=element_blank(),
        axis.ticks.y=element_blank()
  )+
  #各个组织的enhancer数目：
  annotate("text",x=1,y=1.05,label=11137,colour="black",size=rel(3),hjust=0.45)+
  annotate("text",x=2,y=1.05,label=7231,colour="black",size=rel(3),hjust=0.45)+
  annotate("text",x=3,y=1.05,label=1374,colour="black",size=rel(3),hjust=0.45)+
  annotate("text",x=4,y=1.05,label=2727,colour="black",size=rel(3),hjust=0.45)+
  annotate("text",x=5,y=1.05,label=14193,colour="black",size=rel(3),hjust=0.45)+
  annotate("text",x=6,y=1.05,label=6476,colour="black",size=rel(3),hjust=0.45)+
  annotate("text",x=7,y=1.05,label=6485,colour="black",size=rel(3),hjust=0.45)+
  annotate("text",x=8,y=1.05,label=12076,colour="black",size=rel(3),hjust=0.45)+
  annotate("text",x=9,y=1.05,label=5210,colour="black",size=rel(3),hjust=0.45)+
  annotate("text",x=10,y=1.05,label=11256,colour="black",size=rel(3),hjust=0.45)
  
 t(t(apply(table(enh_all_count_spe_melt$tissue,enh_all_count_spe_melt$is_spe),1,sum)))

dev.off()

# 
# 
# 
# ################################################
# #使用组蛋白信号:有副值，没法计算。
# rm(list=ls())
# library(ggplot2)
# setwd("/home/ding/all_cmd/script")
# histone_matrix<-read.table("./enh_statistics/enh_all_cluster_histone_sig.matrix",header = T,stringsAsFactors = F)
# 
# tissues<-names(histone_matrix)[c(-1:-4)]
# #保留两位小数：
# rowsum_sig<-as.matrix(round(rowSums(histone_matrix[,tissues,drop=F]),2),ncol=1)
# #和为0的赋值为0.01
# rowsum_sig[which(rowsum_sig==0)]<-0.01
# for(i in tissues){
#   histone_matrix[,i]<-histone_matrix[,i,drop=F]/rowsum_sig
# }
# 
# #计算香农熵Hs作为组织特异性分值:
# for(i in tissues){
#   histone_matrix[histone_matrix[,i,drop=F]<=0,i]<-log(histone_matrix[histone_matrix[,i,drop=F]<=0,i,drop=F]^histone_matrix[histone_matrix[,i,drop=F]<=0,i,drop=F],2)
#   histone_matrix[histone_matrix[,i,drop=F]>0,i]<-histone_matrix[histone_matrix[,i,drop=F]>0,i,drop=F]*log(histone_matrix[histone_matrix[,i,drop=F]>0,i,drop=F],2)
# 
# }
# 
# ts_matrix<-subset(histone_matrix,select = c("chr","chr_start","chr_end","enhancer","type"))
# ts_matrix$Hs<-(-1)*rowSums(histone_matrix[,tissues,drop=F])
# 


