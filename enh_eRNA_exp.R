#_____________________________________________________________
rm(list=ls())
setwd("/home/ding/all_cmd/script")
library(ggplot2)
# library(pheatmap)
library(reshape2)
#做各个组织bid/unbid eRNA的表达量小提亲图
enh_bid_exp<-read.table("./enh_statistics/enh_bid_exp",header = T,stringsAsFactors = F)
enh_unbid_exp<-read.table("./enh_statistics/enh_unbid_exp",header=T,stringsAsFactors = F)


rownames(enh_bid_exp)<-enh_bid_exp$enhancer
enh_bid_exp<-enh_bid_exp[,names(enh_bid_exp)[5:14]]
enh_bid_exp$type<-"bid"

rownames(enh_unbid_exp)<-enh_unbid_exp$enhancer
enh_unbid_exp<-enh_unbid_exp[,names(enh_unbid_exp)[5:14]]
enh_unbid_exp$type<-"unbid"

enh_all_exp<-rbind(enh_bid_exp,enh_unbid_exp)
enh_all_exp_melt<-melt(enh_all_exp,id.vars = "type",variable.name = "tissue",value.name = "nor_TPM")
names(enh_all_exp_melt)<-c("type","tissue","nor_TPM")
#选择表达量大于0的eRNA做图:
enh_all_exp_melt<-enh_all_exp_melt[enh_all_exp_melt$nor_TPM>1,]
#进行归一化:
 # enh_all_exp_melt$nor_TPM<-enh_all_exp_melt$nor_TPM+1

#分别去除每个组织的out
# nor_TPM_out<-boxplot.stats(enh_all_exp_melt$nor_TPM)$out
# nor_TPM_in<-setdiff(enh_all_exp_melt$nor_TPM,nor_TPM_out)
# enh_all_exp_melt<- enh_all_exp_melt[which(enh_all_exp_melt$nor_TPM %in%  nor_TPM_in),]

a<-data.frame()
#做wilcox.test()检验
for(i in unique(enh_all_exp_melt$tissue)){
  exp_bid_no0<-enh_all_exp_melt[(enh_all_exp_melt$type=="bid")&(enh_all_exp_melt$tissue==i),]$nor_TPM
  exp_unbid_no0<-enh_all_exp_melt[(enh_all_exp_melt$type=="unbid")&(enh_all_exp_melt$tissue==i),]$nor_TPM
  p.value<-as.numeric(wilcox.test(exp_bid_no0,exp_unbid_no0,alternative = "two.sided")$p.value)
  print(paste(i,p.value,sep="  "))
  a_row<-data.frame(tissue=i,p.value)
  a<-rbind(a,a_row)
}

labels_string<-c(0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000)
# ggplot(enh_all_exp_melt,aes(factor(enh_all_exp_melt$tissue),nor_TPM))+geom_violin(aes(fill=type),scale="width")+ylim(0,5)
ggplot(enh_all_exp_melt,aes(factor(enh_all_exp_melt$tissue),nor_TPM))+geom_boxplot(aes(fill=type),show.legend = T)+scale_y_discrete("TPM",breaks=labels_string,labels=labels_string,limits=range(labels_string),expand =c(0.1,0.1) )
ggplot(enh_all_exp_melt,aes(factor(enh_all_exp_melt$tissue),nor_TPM))+geom_boxplot(aes(fill=type),show.legend = T)+scale_y_log10("TPM",breaks=c(0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000),labels=c(0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000))

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_eRNA_exp.png"
CairoPNG(png_path, width = par('din')[1], height = par('din')[2], units='in', dpi=600)

ggplot(enh_all_exp_melt,aes(factor(enh_all_exp_melt$tissue),nor_TPM))+
  geom_boxplot(aes(fill=type),show.legend = T,outlier.size = 0.3)+
  scale_y_log10(breaks=10^(-3:4),labels=paste(breaks=10^(-3:4)))+
  # xlab("tissue")+
  ylab("Expression of eRNAs(TPM)")+
  scale_fill_discrete(limits=c("bid","unbid"),labels=c("2D-eRNA  ","1D-eRNA  "))+
  theme_bw()+
  theme(axis.text=element_text(size=rel(1.1)),
        axis.text.x=element_text(angle = 30,vjust=1,hjust=1),
        axis.text.y=element_text(size=rel(1.1)),
        axis.title.x=element_blank(),
        legend.text=element_text(size=rel(1.1)),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.margin=unit(0,"mm"),
        #legend.position=c(1,1),legend.justification=c(1,1),
        legend.position="top",
        legend.box.spacing = unit(1,"mm")
  )

dev.off()

# rm(list=ls())
# setwd("/home/ding/all_cmd/script")
# library(ggplot2)
# library(pheatmap)
# #做各个组织bid/unbid eRNA的表达量热图
# enh_bid_exp<-read.table("./enh_statistics/enh_bid_exp",header = T,stringsAsFactors = F)
# enh_unbid_exp<-read.table("./enh_statistics/enh_unbid_exp",header=T,stringsAsFactors = F)
# 
# rownames(enh_bid_exp)<-enh_bid_exp$enhancer
# enh_bid_exp<-as.matrix(enh_bid_exp[,names(enh_bid_exp)[5:14]])
# rownames(enh_unbid_exp)<-enh_unbid_exp$enhancer
# enh_unbid_exp<-as.matrix(enh_unbid_exp[,names(enh_unbid_exp)[5:14]])
# 
# enh_all_exp<-rbind(enh_bid_exp,enh_unbid_exp)
# enh_all_exp<-log(enh_all_exp+1,10)
# enh_all_exp[is.infinite(enh_all_exp)]=0
# 
# annotation_col = data.frame(exp_TPM = factor(rep("log(exp_TPM+1)",10)))
# rownames(annotation_col) = colnames(enh_bid_exp)
# 
# annotation_row = data.frame(eRNA = factor(rep(c("bid","unbid"), c(nrow(enh_bid_exp),nrow(enh_unbid_exp)))))
# rownames(annotation_row) = rownames(enh_all_exp)
# 
# ann_colors = list(
#   eRNA = c(bid = "#7570B3", unbid = "#E7298A")
# )
# 
# png_path<-paste("/home/ding/all_cmd/script","/result.mean_signal.pdf",sep="")
# pdf(png_path)
# pheatmap(enh_all_exp,
#          display_numbers=FALSE, 
#          annotation_col = annotation_col, 
#          annotation_row = annotation_row, 
#          annotation_names_row = FALSE,
#          annotation_names_col=TRUE,
#          cluster_row = FALSE,
#          cluster_cols = TRUE,
#          color = colorRampPalette(c("white","red"))(10000),
#          annotation_colors = ann_colors,
#          border_color = "red",
#          scale="none")
# dev.off()

# source("http://bioconductor.org/biocLite.R")
# biocLite("heatmap3")




#################test
# set.seed(3147)
# r<-rnorm(1000)
# summary(r)
# boxplot(r)
# out1<-boxplot.stats(r)
# out<-boxplot.stats(r)$out
# boxplot(r)
# plot(r)
# ok.point <- matrix(c(which(r %in% setdiff(r,out)),setdiff(r,out)),ncol=2)
# out.point <- matrix(c(which(r %in% out), r[which(r %in% out)]), ncol = 2)
# points(out.point,col="red",pch=11)
# points(ok.point,col="blue",pch=16)

