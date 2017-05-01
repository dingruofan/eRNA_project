rm(list=ls())
library(reshape2)
library(ggplot2)
library(stringr)

setwd("/home/ding/all_cmd/script/enh_statistics")
eRNA_anno<-read.table("eRNA_annotate",header = F,stringsAsFactors = F)
eRNA_anno_part<-eRNA_anno[,c(4,10,11)]
names(eRNA_anno_part)<-c("eRNA","gene_name","type")

#lncRNA包括以下几种（根据GENCODE):
# 3prime_overlapping_ncrna
# antisense
# lincRNA
# processed_transcript
# sense_intronic
# sense_overlapping

#将type中的各种类型长非编码类型合并为lncRNA:
eRNA_anno_part[str_detect(eRNA_anno_part$type,"3prime_overlapping_ncrna|lincRNA|antisense|processed_transcript|sense_intronic|sense_overlapping"),"type"]<-"lncRNA"
#观察可知，几乎都是lncRNA和pseudogene。
a<-as.data.frame(t(table(eRNA_anno_part$type)))[,c("Var2","Freq")]
names(a)<-c("gene_type","Freq")
a_other<-data.frame(gene_type="other",Freq=sum(a$Freq)-a[a$gene_type=="lncRNA","Freq"]-a[a$gene_type=="pseudogene","Freq"]-a[a$gene_type=="miRNA","Freq"]-a[a$gene_type=="snRNA","Freq"]-a[a$gene_type=="misc_RNA","Freq"])

eRNA_anno_sta<-rbind(a[a$gene_type=="lncRNA",],a[a$gene_type=="pseudogene",],a_other,a[a$gene_type=="miRNA",],a[a$gene_type=="snRNA",],a[a$gene_type=="misc_RNA",])
eRNA_anno_sta$percent<-paste(round(eRNA_anno_sta$Freq/sum(eRNA_anno_sta$Freq)*100,2),"%",sep="")

#绘制饼图:
eRNA_lables<-paste(eRNA_anno_sta$gene_type,"  (",eRNA_anno_sta$percent,")",sep="")

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/eRNA_annotate_pie.png"
CairoPNG(png_path, width = par('din')[1], height = par('din')[2], units='in', dpi=600)

ggplot(eRNA_anno_sta, aes(x = "" ,y = Freq, fill = gene_type)) +  
  geom_bar(stat = "identity",width = 3)+
  coord_polar(theta = "y")+ labs(x = "", y = "", title = "") +
  # theme(legend.title=element_blank())+
  scale_fill_discrete(breaks=eRNA_anno_sta$gene_type,labels=eRNA_lables)+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        legend.text=element_text(size=rel(1.1)),
        legend.title=element_blank())

dev.off()

a$percent<-paste(round(a$Freq/sum(a$Freq)*100,2),"%",sep="")
a_lables<-paste(a$gene_type,"(",a$percent,")",sep="")
ggplot(a, aes(x = "" ,y = Freq, fill = gene_type)) +   
  geom_bar(stat = "identity",width = 3)+
  coord_polar(theta = "y")+ labs(x = "", y = "", title = "") +
  # theme(legend.title=element_blank())+
  scale_fill_discrete(breaks=a$gene_type,labels=a_lables)+
  theme(axis.text.x=element_blank())
