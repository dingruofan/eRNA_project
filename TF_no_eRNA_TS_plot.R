#从shell读入工作目录
rm(list=ls())
setwd("/home/ding/all_cmd/script")
library(ggplot2)
library(reshape)
library(pheatmap)
library(Cairo)
library(scales)
library(randomcoloR)

TF_no_eRNA_count_exp<-read.table("/home/ding/all_cmd/script/enh_statistics/TF_no_eRNA_count_exp",stringsAsFactors = F,header = T)
TF_bid_count_exp<-read.table("/home/ding/all_cmd/script/enh_statistics/TF_bid_count_exp",stringsAsFactors = F,header = T)
TF_unbid_count_exp<-read.table("/home/ding/all_cmd/script/enh_statistics/TF_unbid_count_exp",stringsAsFactors = F,header = T)
TF_TS_score<-merge(merge(TF_no_eRNA_count_exp,TF_bid_count_exp,by=names(TF_no_eRNA_count_exp),all=T),TF_unbid_count_exp,by=names(TF_no_eRNA_count_exp),all=T)
rm(TF_no_eRNA_count_exp,TF_unbid_count_exp,TF_bid_count_exp)


TF_TS_score<-TF_TS_score[order(TF_TS_score$TS_score),]
TF_TS_score_matrix<-TF_TS_score[,c(2:11)]
rownames(TF_TS_score_matrix)<-TF_TS_score$TF_all_tissue

# ggplot(TF_TS_score,aes(TS_score))+geom_density()
# ggplot(TF_TS_score,aes(x=1,TS_score))+geom_boxplot()

#输出所有转录因子的表达量及TS_score:
write.table(TF_TS_score,file ="./enh_statistics/TF_exp_TSPV",quote = F,sep="\t",col.names = T,row.names = F)

####################################
####################################
#做5个最特异(TSPV最小)的转录因子的各个组织表达量柱状图：

TF_mostTS<-subset(TF_TS_score,TF_all_tissue=="HNF4A"|TF_all_tissue=="FOXA2"|TF_all_tissue=="TFAP2A"|TF_all_tissue=="FOXA1"|TF_all_tissue=="HNF4G"|TF_all_tissue=="TFAP2C")
TF_mostTS_melt<-melt(TF_mostTS,id.vars = c("TF_all_tissue","TS_score"))
names(TF_mostTS_melt)<-c("TF_all_tissue","TS_score","tissue","TF_exp")


png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/TF_exp_mostTS5_filp.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = par('din')[1], height = par('din')[2], units='in', dpi=600)

#如果表达量大于1000，则等于1000:
TF_mostTS_melt1<-TF_mostTS_melt
# TF_mostTS_melt1[TF_mostTS_melt1$TF_exp>4000,"TF_exp"]<-10000

# TF_mostTS_melt1$mask<-0
# TF_mostTS_melt1[TF_mostTS_melt1$TF_exp>4000,"mask"]<-1

# ggplot(TF_mostTS_melt1,aes(x=reorder(TF_all_tissue,TF_exp),y=log(TF_exp+1,2)))+
ggplot(TF_mostTS_melt1,aes(x=reorder(TF_all_tissue,TF_exp),y=(TF_exp+1)))+
  # facet_grid(mask~.,scales="free",space = "free")+
  geom_bar(stat="identity",aes(fill=reorder(tissue,TF_exp)),position = "dodge",width=0.6)+
  # coord_trans(y="log2")+
  # scale_fill_brewer(palette = "Spectral")+
  scale_fill_manual(values = randomColor(10))+
  scale_y_log10(limits=c(10^0,10^5),breaks=c(10^(0:5)),labels=trans_format("log10",math_format(10^.x)),position="left")+
  scale_x_discrete(position="bottom")+
  # expand_limits(y=10^-3)+
  ylab("Expression of TFs")+
  # ylim(0,9100)+
  # coord_cartesian(ylim = c(10^-3, 10^5))+
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(1.3)),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=rel(1.3)),
        legend.margin=margin(0,0,0,0,"mm"),
        legend.position="top",
        legend.text=element_text(size=rel(1.3)),
        legend.box.spacing = unit(1,"mm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.6,"cm")
  )

dev.off()


png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/TF_exp_mostTS5.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = par('din')[1], height = par('din')[2], units='in', dpi=600)


ggplot(TF_mostTS_melt,aes(x=reorder(tissue,TF_exp),y=log(TF_exp+1,2)))+
  geom_bar(stat="identity",aes(fill=reorder(TF_all_tissue,TF_exp)),position = "dodge",width=0.6)+
  scale_fill_brewer(palette = "Spectral")+
  # scale_y_log10(limits=c(10^0,10^5),breaks=c(10^(0:5)),labels=trans_format("log10",math_format(10^.x)),position="left")+
  scale_x_discrete(position="bottom")+
  # expand_limits(y=10^-3)+
  ylab("log Expression of TFs")+
  # coord_cartesian(ylim = c(10^-3, 10^5))+
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(1.1),vjust=1,hjust=1,angle=30),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=rel(1.1)),
        legend.margin=margin(0,0,0,0,"mm"),
        legend.position="top",
        legend.text=element_text(size=rel(1.0)),
        legend.box.spacing = unit(1,"mm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.6,"cm")
  )

dev.off()

# write.table(TF_mostTS,file="/home/ding/桌面/share/TF_mostTS",col.names = T,row.names = T,quote = F,sep="\t")

# TF_TSVT<-read.table("/home/ding/all_cmd/script/enh_statistics/TF_TSVT",header = T,stringsAsFactors = F)
# 
# TF_TS_score[which(TF_TS_score$TF_all_tissue=="HNF4G"),c(2:12)]
# TF_TSVT[which(TF_TSVT$TF_all_tissue=="HNF4G"),c(2:12)]
# 
# 
# max(TF_TS_score[which(TF_TS_score$TF_all_tissue=="HNF4G"),c(2:11)])
# max(TF_TSVT[which(TF_TSVT$TF_all_tissue=="HNF4G"),c(2:11)])

#################################

####################################
####################################

#输出所有转录因子的名称:
write.table(data.frame(TF_all_126=sort(TF_TS_score[,"TF_all_tissue"])),file ="/media/ding/000B49000006264C/eRNA_project/figure/table/database/TF_all_126",quote = F,sep="\t",col.names = T,row.names = F)


#选择1/4分位数作为组织特异性阈值，凡是小于阈值的作为特异性的TF
specific.TF.score<-quantile(TF_TS_score$TS_score)[2]
TF_spe<-TF_TS_score[TF_TS_score$TS_score<=specific.TF.score,]

#生成图片：
png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/TSPV_TF_density.png"
CairoPNG(png_path, width = 6.3, height = 6.3, units='in', dpi=600)

binsize<-diff(range(TF_TS_score$TS_score))/30
ggplot(TF_TS_score,aes(x=TS_score,y=..density..))+
  geom_histogram(binwidth = binsize,alpha=0.3)+
  geom_density(color="red")+
  # geom_vline(xintercept=specific.TF.score,colour="#AAA000", linetype="dotted")+
  scale_x_continuous(breaks=round(as.vector(sort(c(seq(-100,0,10),specific.TF.score))),0))+
  xlab("Tissue-specific values(TSPV) of TFs")+
  ylab("Density")+
  theme_bw()+
  theme(
        axis.text=element_text(size=rel(1.3)),
        axis.text.x=element_text(vjust=1.1),
        axis.title=element_text(size=rel(1.3)),
        legend.text=element_text(size=rel(1.3)),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        #legend.position=c(1,1),legend.justification=c(1,1),
        legend.position="top"
  )
dev.off()

#写入文件:
nrow(TF_spe)
write.table(TF_spe,file ="./enh_statistics/TF_spe",quote = F,sep="\t",col.names = T,row.names = F)


###############################
#按照组织特异性分值从高到低排序后做TF表达量热图，发现并没有什么规律
# annotation_col = data.frame(exp= factor(rep("TF exp",10)))
# rownames(annotation_col) = colnames(TF_TS_score_matrix)
# 
# annotation_row = data.frame(TF=factor(rep("TF",nrow(TF_TS_score_matrix))))
# rownames(annotation_row) = rownames(TF_TS_score_matrix)
# 

# pheatmap(TF_TS_score_matrix,
#          display_numbers=FALSE, 
#          annotation_col = annotation_col, 
#          annotation_row = annotation_row, 
#          annotation_names_row = FALSE,
#          annotation_names_col=TRUE,
#          cluster_row = FALSE,
#          cluster_cols = TRUE,
#          color = colorRampPalette(c("white","red"))(10000),
#          border_color = "red",
#          scale="none")

###############################
#按照组织特异性分值从高到低排序后做T条行图
# rm(list=ls())
# setwd("/home/ding/all_cmd/script")
# library(ggplot2)
# library(reshape)
# library(pheatmap)
# 
# TF_TS_score<-read.table("./enh_statistics/TF_no_eRNA_count_exp",header = T,stringsAsFactors = F)
# 
# 
# TF_TS<-subset(TF_TS_score,select=c("TF_all_tissue","TS_score"))
# #按照组织特异性分值从高到低排序
# TF_TS_order<-TF_TS$TF_all_tissue[order(TF_TS$TS_score)]
# TF_TS$TF_all_tissue<-factor(TF_TS$TF_all_tissue,levels = TF_TS_order)
# 
# ggplot(TF_TS,aes(x=TF_all_tissue,y=TS_score))+geom_bar(stat = "identity")+coord_flip()
# 
# 
# 获得所有与增强子有关的基因:从enh_exp_targetgene_exp
setwd("/home/ding/all_cmd/script/enh_statistics")
enh_target<-read.table("./enh_exp_targetgene_exp",header = F,stringsAsFactors = F)
enh_target_gene=(unique(data.frame(enh_target_gene=sort(enh_target$V10))))
tissues<-unique(data.frame(tissue=enh_target$V13))

write.table(enh_target_gene,"/media/ding/000B49000006264C/eRNA_project/figure/table/database/enh_target_gene",col.names = T,row.names = F,quote = F,append = F)
write.table(tissues,"/media/ding/000B49000006264C/eRNA_project/figure/table/database/tissues",col.names = T,row.names = F,quote = F,append = F)







