#从shell读入工作目录
rm(list=ls())
library(ggplot2)
setwd("/home/ding/all_cmd/script/enh_statistics")

TF_enrich<-read.table("./enh_TF_fold_enrich",header = T,stringsAsFactors = F)
#将enrich为0的换为1,在取log时就成0了：
TF_enrich[TF_enrich[,]==0]<-1

#进行取log,enrich为1的值为0:

TF_enrich$log_eRNA<-log(TF_enrich$eRNA_enrich,10)
TF_enrich$log_no_eRNA<-log(TF_enrich$no_eRNA_enrich,10)

#对颜色进行标记:
TF_enrich$point_color=0
#重点标注以下转录因子:
ann_point<-TF_enrich[TF_enrich$log_no_eRNA==0|TF_enrich$log_eRNA==0,]$TF
TF_enrich[TF_enrich$TF %in% ann_point,"point_color"]<-1

cols<-c("0"="black","1"="red")

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/TF_enrich_isbid.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 7.68, height = 7.1, units='in', dpi=600)

# ggplot(TF_enrich,aes(x=log_eRNA,y=log_no_eRNA))+geom_point()+geom_abline()

P1<-ggplot(TF_enrich,aes(x=log_no_eRNA,y=log_eRNA))+
  geom_point(aes(colour=factor(point_color),size=point_color))+
  scale_color_manual(values=cols)+
  geom_smooth(method="lm",se=F)+
  geom_abline()+
  # xlim(0,7.5)+ylim(0,7.5)+
  scale_size_continuous(range=c(2,3))+
  xlab(expression(paste("log enrichment of ",Enh["noeRNA"]," TFs",sep="")))+
  ylab(expression(paste("log enrichment of ",Enh["eRNA"]," TFs",sep="")))+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1.1)),
        axis.text.x=element_text(size=rel(1.1)),
        legend.position="none"
  )
P1+
  # geom_abline(intercept = -1, slope = 1,color="red")+
  # geom_text(aes(label=TF),size=rel(3),vjust=0.8)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="ESRRA","log_no_eRNA"]+0.1,y=TF_enrich[TF_enrich$TF=="ESRRA","log_eRNA"],label="ESRRA",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="POU5F1","log_no_eRNA"]+0.1,y=TF_enrich[TF_enrich$TF=="POU5F1","log_eRNA"],label="POU5F1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="GRp20","log_no_eRNA"]+0.47,y=TF_enrich[TF_enrich$TF=="GRp20","log_eRNA"],label="GRp20",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="BRF1","log_no_eRNA"]+0.8,y=TF_enrich[TF_enrich$TF=="BRF1","log_eRNA"],label="BRF1",colour="red",size=rel(2.5),hjust=0)+
  annotate("text",x=TF_enrich[TF_enrich$TF=="WRNIP1","log_no_eRNA"],y=TF_enrich[TF_enrich$TF=="WRNIP1","log_eRNA"]+0.1,label="WRNIP1",colour="red",size=rel(2.5),vjust=0)


  dev.off()
  
#重点标注以下转录因子:
b<-TF_enrich[TF_enrich$log_no_eRNA==0|TF_enrich$log_eRNA==0,]
b<-b[order(b[,"log_eRNA"],decreasing = T),]
b

