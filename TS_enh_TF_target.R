#统计各个组织的TF.
rm(list=ls())
library(ggplot2)
library(plyr)
library(reshape2)
setwd("/home/ding/all_cmd/script/enh_statistics")

#先是TS enh和 TS TF的图:
spe_TF_enh<-read.table("./spe_TF_enh_tissue_2",header = F,stringsAsFactors = F)

#载入TF的特异性score,以及TSVT
TF_TS_score<-read.table("./TF_TSVT",header = T,stringsAsFactors = F)
#对列名重新命名:
names(TF_TS_score)[1]<-"TF"
names(spe_TF_enh)<-c("enh","type","TF","tissue")
#处理数据：
TF_TS_score_melt<-melt(TF_TS_score[,-12],id.vars = "TF",variable_name = "tissue")
names(TF_TS_score_melt)[2:3]<-c("tissue","TSVT")

#_______________________________
#统计富集程度:
spe_TF_enh_table<-data.frame(t(table(spe_TF_enh[,c("type","tissue")])))
#对type排序：
spe_TF_enh_table$type<-factor(spe_TF_enh_table$type,levels=c("bid","unbid","no_eRNA"))
#TS enh和其对应的TS TF的富集程度:
ggplot(spe_TF_enh_table,aes(x=reorder(tissue,Freq,sum),y=Freq,fill=type,order=type))+
  geom_bar(stat = "identity",position = "stack")


#_______________________________
#另外一种点图:
spe_TF_enh_table2<-data.frame(t(table(spe_TF_enh[,c("TF","tissue")])))
#增加TF的特异性score来表示颜色：
spe_TF_enh_table2<-merge(spe_TF_enh_table2,TF_TS_score_melt,by=c("TF","tissue"))
#排序：
spe_TF_enh_table2$TF<-reorder(spe_TF_enh_table2$TF,spe_TF_enh_table2$Freq)
spe_TF_enh_table2$tissue<-reorder(spe_TF_enh_table2$tissue,spe_TF_enh_table2$Freq)
#TF富集程度使用log(2)归一化，TF富集程度小与等于2的全算为2:
spe_TF_enh_table2[spe_TF_enh_table2$Freq<=2&spe_TF_enh_table2$Freq!=0,"Freq"]<-2
#Freq为0的替换为NA:
spe_TF_enh_table2[spe_TF_enh_table2$Freq==0,"TSVT"]<-NA
spe_TF_enh_table2[spe_TF_enh_table2$Freq==0,"Freq"]<-NA


png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/TS_enh_TF_point.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 4.69, height = 7.2, units='in', dpi=600)

ggplot(spe_TF_enh_table2,aes(x=tissue,y=TF))+
  geom_point(na.rm = T,shape=21,stat="identity",aes(fill=TSVT,size=log(Freq,2)))+
  scale_fill_gradient(limits=c(-20,0),low="white",high="red")+
  scale_size(expression(log["2"](TF~enrichment)),breaks=seq(1,10,2),range=c(1,6))+
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(1.1),vjust=1,hjust=1,angle=30),
        axis.title=element_blank(),
        legend.text=element_text(size=rel(1.1)),
        panel.grid.major = element_blank()
  )

dev.off()

#临时输出富集程度:#############
spe_TF_enh_table2$enrich<-log(spe_TF_enh_table2$Freq,2)
spe_TF_enh_table2[which(is.na(spe_TF_enh_table2$enrich)),"enrich"]<-0

TF_select="TFAP2A"
spe_TF_enh_table2[which(spe_TF_enh_table2$TF==TF_select),,drop=F][order(spe_TF_enh_table2[which(spe_TF_enh_table2$TF==TF_select),,drop=F]$enrich),]
quantile(spe_TF_enh_table2[which(spe_TF_enh_table2$TF==TF_select),]$enrich,probs = seq(0,1,0.2))

#############

#临时输出每个组织最特异的转录因子:#############


tissues<-unique(spe_TF_enh_table2$tissue)
a<-data.frame()
for(i in tissues){
  a_tmp<-spe_TF_enh_table2[which(spe_TF_enh_table2$tissue==i),]
  a_tmp<-a_tmp[complete.cases(a_tmp),]
  a_tmp<-a_tmp[order(a_tmp$TSVT,decreasing = T),][1,]
  a<-rbind(a,a_tmp)
}
write.table(a,"/home/ding/all_cmd/script/enh_statistics/10_mostTS_TF",row.names = F,col.names = T,quote = F)

#然后获得与这10个转录因子相关的TS enhancer/TS target gene:
TS_eTG<-read.table("/media/ding/000B49000006264C/eRNA_project/figure/table/database/TS_enh_TF_targetgene_3",header = T,stringsAsFactors = F)
TS_eTG_a<-merge(a,TS_eTG,by=c("tissue","TF"))
write.table(TS_eTG_a,"/media/ding/000B49000006264C/eRNA_project/figure/table/database/TS_enh_TF_targetgene_10_mostTS_TF",row.names = F,col.names = T,quote = F)

#每个组织有结合位点的最特异的转录因子:
tissues<-unique(TS_eTG$tissue)
a<-data.frame()
for(i in tissues){
  a_tmp<-TS_eTG[TS_eTG$tissue==i,]
  a_tmp_TF<-a_tmp[order(a_tmp$TF_TS_score),"TF"][1]
  a_tmp_2<-a_tmp[a_tmp$tissue==i&a_tmp$TF==a_tmp_TF,]
  a<-rbind(a,a_tmp_2)
}

#############

#只做6个转录因子在10个组织中的表达量的柱状图：
TF_mostTS<-subset(spe_TF_enh_table2,TF=="FOXA1"|TF=="FOXA2"|TF=="HNF4A"|TF=="HNF4G"|TF=="TFAP2A"|TF=="TFAP2C")
TF_mostTS$nor_Freq<-log(TF_mostTS$Freq,2)

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/TS_enh_TF_point_mostTS5.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = par('din')[1], height = par('din')[2], units='in', dpi=600)

ggplot(TF_mostTS,aes(x=tissue,y=nor_Freq))+
  geom_bar(stat="identity",aes(fill=TF),position = "dodge",width=0.6)+
  scale_fill_brewer(palette = "Dark2")+
  ylab("TF enrichment")+
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(1.1),vjust=1,hjust=1,angle=30),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=rel(1.1)),
        legend.margin=margin(0,0,0,0,"mm"),
        legend.position="top",
        legend.text=element_text(size=rel(1.1)),
        legend.box.spacing = unit(1,"mm"),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.box="horizontal",
        legend.key.size=unit(0.6,"cm")
  )

dev.off()






#_______________________________
#_______________________________
#_______________________________
#然后是TS enh和 TS target gene的图:
rm(list=ls())
library(ggplot2)
library(plyr)
library(reshape2)
setwd("/home/ding/all_cmd/script/enh_statistics")

spe_enh_target<-read.table("./TS_enh_targetgene",header = F,stringsAsFactors = F)
names(spe_enh_target)<-c("enh","gene_id","gene_name","gene_exp","type","tissue","DPM","isspe")


#_______________________________
#统计targetgene富集程度:
spe_enh_target_table<-data.frame(t(table(spe_enh_target[!duplicated(spe_enh_target[,c("gene_name","tissue")]),c("tissue")])))
names(spe_enh_target_table)[2]<-"tissue"
#TS enh和其对应的TS gene的富集程度:
ggplot(spe_enh_target_table,aes(x=reorder(tissue,Freq),y=Freq))+geom_bar(stat = "identity",position = "stack")+theme_linedraw()

#做每个组织特异性enh调控的特异性gene的表达量箱线图：
spe_enh_target_exp<-spe_enh_target[!duplicated(spe_enh_target[,c("gene_name","gene_exp","tissue")]),c("gene_name","gene_exp","tissue")]
ggplot(spe_enh_target_exp,aes(x=reorder(tissue,gene_exp,median),y=gene_exp))+scale_y_log10(breaks=10^(-1:4),labels=10^(-1:4))+geom_boxplot()

#做点图:
spe_enh_target_1<-spe_enh_target[!duplicated(spe_enh_target[,c("gene_name","gene_exp","tissue","DPM")]),c("gene_name","gene_exp","tissue","DPM")]
ggplot(spe_enh_target_1,aes(y=gene_name,x=tissue))+geom_point(stat="identity")


#########################################
#_______________________________
#特异性的TF调控的TS gene和TS TF:
rm(list=ls())
library(ggplot2)
library(plyr)
library(reshape2)
TS_enh_TF_target<-read.table("./TS_enh_TF_targetgene_2",header = F,stringsAsFactors = F)
names(TS_enh_TF_target)<-c("enh","type","TF","tissue","gene_name","RPKM","DPM")

#有以下TF循环调控靶gene:
TS_enh_TF_target[TS_enh_TF_target$gene_name==TS_enh_TF_target$TF,]
#图列为TF，行为组织，点表示被调控gene的表达量和特异性
TS_enh_TF_target_1<-unique(TS_enh_TF_target[,c("TF","tissue","gene_name","RPKM","DPM")])

#所有特异性gene做图，
tissue_order<-as.data.frame(t(table(TS_enh_TF_target_1$tissue)))[,c(2,3)]
names(tissue_order)<-c("tissue","Freq")
TS_enh_TF_target_2<-merge(TS_enh_TF_target_1,tissue_order,by="tissue")
ggplot(TS_enh_TF_target_2,aes(x=reorder(tissue,Freq),y=reorder(gene_name,DPM),fill=log(RPKM+1,10)))+geom_raster()+scale_fill_gradient2(low="white",high="red")


#统计每个组织中每个转录因子调控的gene的数目:
rm(list=ls())
library(ggplot2)
library(plyr)
library(reshape2)
TS_enh_TF_target_3<-read.table("./TS_enh_TF_targetgene_3",header = F,stringsAsFactors = F)
names(TS_enh_TF_target_3)<-c("enh","type","TF","tissue","gene_name","RPKM","DPM","TSPV")

TS_enh_TF_target_3_tmp<-TS_enh_TF_target_3[!duplicated(TS_enh_TF_target_3[,c("TF","tissue","TSPV")]),c("TF","tissue","TSPV")]
TS_enh_TF_target_3_tmp$gene_num<-0
for(i in unique(TS_enh_TF_target_3$tissue)){
  for(j in unique(TS_enh_TF_target_3$TF)){
    if(length(unique(TS_enh_TF_target_3[TS_enh_TF_target_3$tissue==i&TS_enh_TF_target_3$TF==j,"gene_name"]))>0){
    TS_enh_TF_target_3_tmp[TS_enh_TF_target_3_tmp$tissue==i&TS_enh_TF_target_3_tmp$TF==j,"gene_num"] <- length(unique(TS_enh_TF_target_3[TS_enh_TF_target_3$tissue==i&TS_enh_TF_target_3$TF==j,"gene_name"]))
   }
  }
}

#统计每个组织TS enh/TF/gene的数目:

ggplot(TS_enh_TF_target_3_tmp,aes(x=reorder(tissue,gene_num),y=reorder(TF,TSPV)))+geom_point(aes(size=gene_num))+scale_size_continuous(breaks=c(1,5,10,20,50,80,100))

######################
for(i in unique(TS_enh_TF_target_3_tmp$tissue)){
  TS_enh_TF_target_3_tmp[TS_enh_TF_target_3_tmp$tissue==i,"TF_num"]<-length(unique(TS_enh_TF_target_3[TS_enh_TF_target_3$tissue==i,"TF"]))
  TS_enh_TF_target_3_tmp[TS_enh_TF_target_3_tmp$tissue==i,"enh_num"]<-length(unique(TS_enh_TF_target_3[TS_enh_TF_target_3$tissue==i,"enh"]))
  TS_enh_TF_target_3_tmp[TS_enh_TF_target_3_tmp$tissue==i,"gene_num"]<-length(unique(TS_enh_TF_target_3[TS_enh_TF_target_3$tissue==i,"gene_name"]))
}
TS_enh_TF_target_4<-TS_enh_TF_target_3_tmp[!duplicated(TS_enh_TF_target_3_tmp[,c("tissue","gene_num","TF_num","enh_num")]),c("tissue","gene_num","TF_num","enh_num")]

TS_enh_TF_target_4_melt<-melt(TS_enh_TF_target_4,id="tissue")
names(TS_enh_TF_target_4_melt)<-c("tissue","type","Freq")
TS_enh_TF_target_4_melt$type<-factor(TS_enh_TF_target_4_melt$type,levels=c("TF_num","gene_num","enh_num"))
ggplot(TS_enh_TF_target_4_melt,aes(x=reorder(tissue,Freq,median),y=Freq))+geom_bar(aes(fill=type),stat="identity",position = "dodge")








