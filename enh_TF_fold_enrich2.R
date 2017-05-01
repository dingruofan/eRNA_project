#从shell读入工作目录
rm(list=ls())
library(ggplot2)
setwd("/home/ding/all_cmd/script/enh_statistics")

TF_enrich<-read.table("./enh_TF_fold_enrich2",header = T,stringsAsFactors = F)
#将enrich为0的换为1,在取log时就成0了：
TF_enrich[TF_enrich[,]==0]<-1

#进行取log,enrich为1的值为0:

TF_enrich$log_spe<-log(TF_enrich$spe_enrich,10)
TF_enrich$log_uni<-log(TF_enrich$uni_enrich,10)

#对颜色进行标记:
TF_enrich$point_color=0
#重点标注以下转录因子:
ann_point<-TF_enrich[TF_enrich$log_spe==0|TF_enrich$uni_enrich==0,]$TF
TF_enrich[TF_enrich$TF %in% ann_point,"point_color"]<-1

cols<-c("0"="black","1"="red")

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/TF_enrich_isspe.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 7.68, height = 7.1, units='in', dpi=600)

# ggplot(TF_enrich,aes(x=log_spe,y=log_uni))+geom_point()+geom_abline()

P1<-ggplot(TF_enrich,aes(x=log_spe ,y=log_uni))+
  geom_point(aes(colour=factor(point_color),size=point_color))+
  scale_color_manual(values=cols)+
  geom_abline()+
  xlim(0,5.5)+ylim(0,5.5)+
  scale_size_continuous(range=c(2,3))+
  xlab(expression(paste("log enrichment of ",Enh["TS"]," TFs",sep="")))+
  ylab(expression(paste("log enrichment of ",Enh["UE"]," TFs",sep="")))+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1.1)),
        axis.text.x=element_text(size=rel(1.1)),
        legend.position="none"
  )


P1+
  # geom_abline(intercept = -1, slope = 1,color="red")+
  # geom_text(aes(label=TF),size=rel(3),vjust=0.8)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="JUN","log_spe"]+0.1,y=TF_enrich[TF_enrich$TF=="JUN","log_uni"],label="JUN",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="POU2F2","log_spe"]+0.1,y=TF_enrich[TF_enrich$TF=="POU2F2","log_uni"],label="POU2F2",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="ZBTB7A","log_spe"]+0.1,y=TF_enrich[TF_enrich$TF=="ZBTB7A","log_uni"],label="ZBTB7A",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="HNF4G","log_spe"]+0.5,y=TF_enrich[TF_enrich$TF=="HNF4G","log_uni"],label="HNF4G",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="HNF4A","log_spe"]+0.1,y=TF_enrich[TF_enrich$TF=="HNF4A","log_uni"],label="HNF4A",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="EZH2","log_spe"]+0.45,y=TF_enrich[TF_enrich$TF=="EZH2","log_uni"],label="EZH2",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="BCL3","log_spe"]+0.83,y=TF_enrich[TF_enrich$TF=="BCL3","log_uni"],label="BCL3",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="PBX3","log_spe"]+1.2,y=TF_enrich[TF_enrich$TF=="PBX3","log_uni"],label="PBX3",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="SIX5","log_spe"]+1.5,y=TF_enrich[TF_enrich$TF=="SIX5","log_uni"],label="SIX5",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="ATF1","log_spe"]+1.8,y=TF_enrich[TF_enrich$TF=="ATF1","log_uni"],label="ATF1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="MYBL2","log_spe"]+2.1,y=TF_enrich[TF_enrich$TF=="MYBL2","log_uni"],label="MYBL2",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="ETS1","log_spe"]+0.9,y=TF_enrich[TF_enrich$TF=="ETS1","log_uni"]+0.08,label="ETS1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="ZBTB33","log_spe"]+0.1,y=TF_enrich[TF_enrich$TF=="ZBTB33","log_uni"]+0.1,label="ZBTB33",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="BRCA1","log_spe"]+0.5,y=TF_enrich[TF_enrich$TF=="BRCA1","log_uni"]+0.15,label="BRCA1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="KAP1","log_spe"]+0.5,y=TF_enrich[TF_enrich$TF=="KAP1","log_uni"]+0.01,label="KAP1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="SIN3AK20","log_spe"]+0.1,y=TF_enrich[TF_enrich$TF=="SIN3AK20","log_uni"]+0.06,label="SIN3AK20",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="FOXM1","log_spe"]+0.1,y=TF_enrich[TF_enrich$TF=="FOXM1","log_uni"],label="FOXM1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="JUNB","log_spe"]+0.5,y=TF_enrich[TF_enrich$TF=="JUNB","log_uni"],label="JUNB",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="RPC155","log_spe"]+0.1,y=TF_enrich[TF_enrich$TF=="RPC155","log_uni"],label="RPC155",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="ZKSCAN1","log_spe"]+0.5,y=TF_enrich[TF_enrich$TF=="ZKSCAN1","log_uni"],label="ZKSCAN1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="BCLAF1","log_spe"]+0.9,y=TF_enrich[TF_enrich$TF=="BCLAF1","log_uni"],label="BCLAF1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="THAP1","log_spe"]+1.3,y=TF_enrich[TF_enrich$TF=="THAP1","log_uni"],label="THAP1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="GTF2F1","log_spe"]+1.6,y=TF_enrich[TF_enrich$TF=="GTF2F1","log_uni"],label="GTF2F1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="BCL11A","log_spe"]+2,y=TF_enrich[TF_enrich$TF=="BCL11A","log_uni"],label="BCL11A",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="UBTF","log_spe"]+2.4,y=TF_enrich[TF_enrich$TF=="UBTF","log_uni"],label="UBTF",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="IRF4","log_spe"]+0.1,y=TF_enrich[TF_enrich$TF=="IRF4","log_uni"],label="IRF4",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="BDP1","log_spe"]+0.4,y=TF_enrich[TF_enrich$TF=="BDP1","log_uni"],label="BDP1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="GTF3C2","log_spe"]+0.8,y=TF_enrich[TF_enrich$TF=="GTF3C2","log_uni"],label="GTF3C2",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="HDAC2","log_spe"]+1.2,y=TF_enrich[TF_enrich$TF=="HDAC2","log_uni"],label="HDAC2",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="PRDM1","log_spe"]+1.5,y=TF_enrich[TF_enrich$TF=="PRDM1","log_uni"],label="PRDM1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="ELK1","log_spe"]+0.4,y=TF_enrich[TF_enrich$TF=="ELK1","log_uni"],label="ELK1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="IKZF1","log_spe"]+0.1,y=TF_enrich[TF_enrich$TF=="IKZF1","log_uni"],label="IKZF1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="STAT2","log_spe"]+0.8,y=TF_enrich[TF_enrich$TF=="STAT2","log_uni"],label="STAT2",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="SP2","log_spe"]+1.3,y=TF_enrich[TF_enrich$TF=="SP2","log_uni"],label="SP2",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="RBBP5","log_spe"]+0.1,y=TF_enrich[TF_enrich$TF=="RBBP5","log_uni"],label="RBBP5",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="NR2C2","log_spe"]+0.5,y=TF_enrich[TF_enrich$TF=="NR2C2","log_uni"],label="NR2C2",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="GTF2B","log_spe"]+0.85,y=TF_enrich[TF_enrich$TF=="GTF2B","log_uni"],label="GTF2B",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="FOXP2","log_spe"]+1.2,y=TF_enrich[TF_enrich$TF=="FOXP2","log_uni"],label="FOXP2",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="CEBPD","log_spe"]+0.1,y=TF_enrich[TF_enrich$TF=="CEBPD","log_uni"],label="CEBPD",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="CTCFL","log_spe"]+0.5,y=TF_enrich[TF_enrich$TF=="CTCFL","log_uni"],label="CTCFL",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="MBD4","log_spe"]+0.85,y=TF_enrich[TF_enrich$TF=="MBD4","log_uni"],label="MBD4",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="MEF2C","log_spe"]+1.2,y=TF_enrich[TF_enrich$TF=="MEF2C","log_uni"],label="MEF2C",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="IRF3","log_spe"]+1.5,y=TF_enrich[TF_enrich$TF=="IRF3","log_uni"],label="IRF3",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="SETDB1","log_spe"]+1.8,y=TF_enrich[TF_enrich$TF=="SETDB1","log_uni"],label="SETDB1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="SIN3A","log_spe"]+2.2,y=TF_enrich[TF_enrich$TF=="SIN3A","log_uni"],label="SIN3A",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="ZNF217","log_spe"]+2.5,y=TF_enrich[TF_enrich$TF=="ZNF217","log_uni"],label="ZNF217",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="SMARCC1","log_spe"]+0.1,y=TF_enrich[TF_enrich$TF=="SMARCC1","log_uni"],label="SMARCC1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="BRF2","log_spe"]+0.6,y=TF_enrich[TF_enrich$TF=="BRF2","log_uni"],label="BRF2",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="ZEB1","log_spe"]+1,y=TF_enrich[TF_enrich$TF=="ZEB1","log_uni"],label="ZEB1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="NFE2","log_spe"]+1.3,y=TF_enrich[TF_enrich$TF=="NFE2","log_uni"],label="NFE2",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="HSF1","log_spe"]+1.6,y=TF_enrich[TF_enrich$TF=="HSF1","log_uni"],label="HSF1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="SP4","log_spe"]+2,y=TF_enrich[TF_enrich$TF=="SP4","log_uni"],label="SP4",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="ZNF274","log_spe"]+2.3,y=TF_enrich[TF_enrich$TF=="ZNF274","log_uni"],label="ZNF274",colour="red",size=rel(2.5),hjust=0)

dev.off()


b<-TF_enrich[TF_enrich$log_spe==0|TF_enrich$log_uni==0,]
b<-b[order(b[,"log_uni"],decreasing = T),]
b

  