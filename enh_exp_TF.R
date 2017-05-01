
rm(list=ls())
library(ggplot2)
library(reshape2)
library(Cairo)
setwd("/home/ding/all_cmd/script/enh_statistics/")
enh_bid_exp_TF<-read.table("./enh_bid_exp_TF_all",header = T,stringsAsFactors = F)
enh_unbid_exp_TF<-read.table("./enh_unbid_exp_TF_all",header = T,stringsAsFactors = F)
enh_no_eRNA_exp_TF<-read.table("./enh_no_eRNA_exp_TF_all",header = T,stringsAsFactors = F)
enh_random_exp_TF<-read.table("./enh_random_exp_TF_all",header = T,stringsAsFactors = F)


names(enh_bid_exp_TF)<-c("enhancer","eRNA_exp_TPM","TFBS_num","tissue")
names(enh_unbid_exp_TF)<-c("enhancer","eRNA_exp_TPM","TFBS_num","tissue")
names(enh_no_eRNA_exp_TF)<-c("enhancer","eRNA_exp_TPM","TFBS_num","tissue")
names(enh_random_exp_TF)<-c("enhancer","eRNA_exp_TPM","TFBS_num","tissue")

enh_all_exp_TF<-rbind(enh_bid_exp_TF,enh_unbid_exp_TF,enh_no_eRNA_exp_TF,enh_random_exp_TF)
enh_all_exp_TF$type<-c(rep("bid",nrow(enh_bid_exp_TF)),rep("unbid",nrow(enh_unbid_exp_TF)),rep("no_eRNA",nrow(enh_no_eRNA_exp_TF)),rep("random",nrow(enh_random_exp_TF)))

# cor.test(enh_bid_exp_TF$eRNA_exp_TPM,enh_bid_exp_TF$TFBS_num,method=c("pearson"))
# cor.test(enh_unbid_exp_TF$eRNA_exp_TPM,enh_unbid_exp_TF$TFBS_num,method=c("pearson"))
enh_all_exp_TF<-enh_all_exp_TF[!duplicated(enh_all_exp_TF[,c(1,3,4,5)]),c(1,3,4,5)]
#!!!！！!!！!注意：使用log(TFBS_num)作为转录因子富集程度:
# enh_all_exp_TF$TFBS_num<-log(enh_all_exp_TF,2)

#排序：
enh_all_exp_TF$type<-factor(enh_all_exp_TF$type,levels=c("bid","unbid","no_eRNA","random"))
#######################################

#bid/unbid/no_eRNA转录因子结合位点数目统计箱线图
# par('din')
png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_isbid_TF_enrich.png"
CairoPNG(png_path, width = 5.5, height =5.16 , units='in', dpi=600)

ggplot(enh_all_exp_TF,aes(x=reorder(tissue,TFBS_num),y=TFBS_num))+
  geom_boxplot(aes(fill=type),show.legend = T,outlier.size = 0,outlier.stroke = 0)+
  ylim(0,40)+
  xlab("tissue")+
  ylab("TF enrichment")+
  scale_fill_discrete(limits=c("bid","unbid","no_eRNA","random"),
                      labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"]),"random"))+
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
  )

dev.off()

enh_TF<-enh_all_exp_TF[enh_all_exp_TF$type!="random","TFBS_num"]
random_TF<-enh_all_exp_TF[enh_all_exp_TF$type=="random","TFBS_num"]
t.test(enh_TF,random_TF,alternative = "greater")


#做检验：
a<-data.frame()
for(i in unique(enh_all_exp_TF$tissue)){
  bid_unbid_p<-(t.test(enh_all_exp_TF[enh_all_exp_TF$tissue==i&enh_all_exp_TF$type=="bid","TFBS_num"],enh_all_exp_TF[enh_all_exp_TF$tissue==i&enh_all_exp_TF$type=="unbid","TFBS_num"],alternative = "greater")$p.value)
  bid_no_eRNA_p<-(t.test(enh_all_exp_TF[enh_all_exp_TF$tissue==i&enh_all_exp_TF$type=="bid","TFBS_num"],enh_all_exp_TF[enh_all_exp_TF$tissue==i&enh_all_exp_TF$type=="no_eRNA","TFBS_num"],alternative = "greater")$p.value)
  bid_random_p<-(t.test(enh_all_exp_TF[enh_all_exp_TF$tissue==i&enh_all_exp_TF$type=="bid","TFBS_num"],enh_all_exp_TF[enh_all_exp_TF$tissue==i&enh_all_exp_TF$type=="random","TFBS_num"],alternative = "greater")$p.value)
  unbid_no_eRNA_p<-(t.test(enh_all_exp_TF[enh_all_exp_TF$tissue==i&enh_all_exp_TF$type=="unbid","TFBS_num"],enh_all_exp_TF[enh_all_exp_TF$tissue==i&enh_all_exp_TF$type=="no_eRNA","TFBS_num"],alternative = "greater")$p.value)
  unbid_random_p<-(t.test(enh_all_exp_TF[enh_all_exp_TF$tissue==i&enh_all_exp_TF$type=="unbid","TFBS_num"],enh_all_exp_TF[enh_all_exp_TF$tissue==i&enh_all_exp_TF$type=="random","TFBS_num"],alternative = "greater")$p.value)
  no_eRNA_random_p<-(t.test(enh_all_exp_TF[enh_all_exp_TF$tissue==i&enh_all_exp_TF$type=="no_eRNA","TFBS_num"],enh_all_exp_TF[enh_all_exp_TF$tissue==i&enh_all_exp_TF$type=="random","TFBS_num"],alternative = "greater")$p.value)
  a_row<-data.frame(tissue=i,Enh2D_Enh1D=bid_unbid_p,Enh2D_EnhnoeRNA=bid_no_eRNA_p,Enh2D_random=bid_random_p,Enh1D_EnhnoeRNA=unbid_no_eRNA_p,Enh1D_random=unbid_random_p,EnhnoRNA_random=no_eRNA_random_p)
  a<-rbind(a,a_row)
  print(paste(i,"  2D-1D:",bid_unbid_p,"  2D-no-eRNA:",bid_no_eRNA_p,"    2D-random:",bid_random_p,"  1D-no_eRNA:",unbid_no_eRNA_p,"  1D-random:",unbid_random_p,"  no-eRNA-random:",no_eRNA_random_p,sep=" "))
}

a<-data.frame()
for(i in unique(enh_all_exp_TF$tissue)){
  enh_random_p<-(t.test(enh_all_exp_TF[enh_all_exp_TF$tissue==i,"TFBS_num"],enh_all_exp_TF[enh_all_exp_TF$tissue==i&enh_all_exp_TF$type=="random","TFBS_num"],alternative = "greater")$p.value)
  a_row<-data.frame(tissue=i,enh_random=enh_random_p)
  a<-rbind(a,a_row)
  print(paste(i,"  enh-random:",enh_random_p))
}


#获得spe的信息，并做图:
enh_all_count_spe<-read.table("/home/ding/all_cmd/script/enh_statistics/enh_all_count_spe",header = T,stringsAsFactors = F)
enh_all_count_spe<-enh_all_count_spe[,c("enhancer","type","is_spe")]
#合并spe列
enh_all_exp_TF_spe<-merge(enh_all_exp_TF,enh_all_count_spe,by=c("enhancer","type"),all.x = T)
#将random行的is_spe改为random:
enh_all_exp_TF_spe[is.na(enh_all_exp_TF_spe$is_spe),"is_spe"]<-"random"

#排序：
enh_all_exp_TF_spe$is_spe<-factor(enh_all_exp_TF_spe$is_spe,levels=c("uni","other","spe","random"))
#######################################

#bid/unbid/no_eRNA转录因子结合位点数目统计箱线图
# par('din')
png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_isspe_TF_enrich.png"
CairoPNG(png_path, width = 5.5, height =5.16 , units='in', dpi=600)

ggplot(enh_all_exp_TF_spe,aes(x=reorder(tissue,TFBS_num),y=TFBS_num))+
  geom_boxplot(aes(fill=is_spe),show.legend = T,outlier.size = 0,outlier.stroke = 0)+
  ylim(0,30)+
  xlab("tissue")+
  ylab("TF enrichment")+
  scale_fill_discrete(limits=c("spe","other","uni","random"),
                      labels=c(expression(Enh["TS"]),expression(Enh["Oth"]),expression(Enh["UE"]),"random"))+
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
  )
dev.off()


#做检验：
a<-data.frame()
for(i in unique(enh_all_exp_TF_spe$tissue)){
  UE_Oth_p<-(t.test(enh_all_exp_TF_spe[enh_all_exp_TF_spe$tissue==i&enh_all_exp_TF_spe$is_spe=="uni","TFBS_num"],enh_all_exp_TF_spe[enh_all_exp_TF_spe$tissue==i&enh_all_exp_TF_spe$is_spe=="other","TFBS_num"],alternative = "greater")$p.value)
  UE_TS_p<-(t.test(enh_all_exp_TF_spe[enh_all_exp_TF_spe$tissue==i&enh_all_exp_TF_spe$is_spe=="uni","TFBS_num"],enh_all_exp_TF_spe[enh_all_exp_TF_spe$tissue==i&enh_all_exp_TF_spe$is_spe=="spe","TFBS_num"],alternative = "greater")$p.value)
  UE_random_p<-(t.test(enh_all_exp_TF_spe[enh_all_exp_TF_spe$tissue==i&enh_all_exp_TF_spe$is_spe=="uni","TFBS_num"],enh_all_exp_TF_spe[enh_all_exp_TF_spe$tissue==i&enh_all_exp_TF_spe$is_spe=="random","TFBS_num"],alternative = "greater")$p.value)
  Oth_TS_p<-(t.test(enh_all_exp_TF_spe[enh_all_exp_TF_spe$tissue==i&enh_all_exp_TF_spe$is_spe=="other","TFBS_num"],enh_all_exp_TF_spe[enh_all_exp_TF_spe$tissue==i&enh_all_exp_TF_spe$is_spe=="spe","TFBS_num"],alternative = "greater")$p.value)
  Oth_random_p<-(t.test(enh_all_exp_TF_spe[enh_all_exp_TF_spe$tissue==i&enh_all_exp_TF_spe$is_spe=="other","TFBS_num"],enh_all_exp_TF_spe[enh_all_exp_TF_spe$tissue==i&enh_all_exp_TF_spe$is_spe=="random","TFBS_num"],alternative = "greater")$p.value)
  TS_random_p<-(t.test(enh_all_exp_TF_spe[enh_all_exp_TF_spe$tissue==i&enh_all_exp_TF_spe$is_spe=="spe","TFBS_num"],enh_all_exp_TF_spe[enh_all_exp_TF_spe$tissue==i&enh_all_exp_TF_spe$is_spe=="random","TFBS_num"],alternative = "greater")$p.value)
  a_row<-data.frame(tissue=i,EnhUE_EnhOth=UE_Oth_p,EnhUE_EnhTS=UE_TS_p,EnhUE_random=UE_random_p,
                    EnhOth_EnhTS=Oth_TS_p,EnhOth_random=Oth_random_p,EnhTS_random=TS_random_p)
  a<-rbind(a,a_row)
  print(paste(i,"  UE-Oth:",UE_Oth_p,"  UE-TS:",UE_TS_p,"    UE-random:",UE_random_p,"  Oth-TS:",Oth_TS_p,"  Oth-random:",Oth_random_p,"  TS-random:",TS_random_p,sep=" "))
}


#比较eRNAs-ENh和no-eRNAs Enh的TF富集程度：
enh_all_exp_TF_2<-enh_all_exp_TF
enh_all_exp_TF_2$type<-as.character(enh_all_exp_TF_2$type)
enh_all_exp_TF_2[enh_all_exp_TF_2$type=="bid"|enh_all_exp_TF_2$type=="unbid","type"]<-"eRNAs-Enh"

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/enh_iseRNAs_TF_enrich.png"
CairoPNG(png_path, width = 5.5, height =5.16 , units='in', dpi=600)

ggplot(enh_all_exp_TF_2,aes(x=reorder(tissue,TFBS_num),y=TFBS_num))+
  geom_boxplot(aes(fill=type),show.legend = T,outlier.size = 0,outlier.stroke = 0)+
  ylim(0,30)+
  xlab("tissue")+
  ylab("TF enrichment")+
  scale_fill_discrete(limits=c("eRNAs-Enh","no_eRNA","random"),labels=c(expression(Enh["eRNA"]),expression(Enh["no-eRNA"]),"random"))+
  theme_bw()+
  theme(
        axis.text.x=element_text(angle = 30,vjust=1,hjust=1),
        axis.title.x=element_blank(),
        legend.text=element_text(size=rel(1.1)),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.margin=unit(0,"mm"),
        legend.box.spacing = unit(0.1,"mm"),
        #legend.position=c(1,1),legend.justification=c(1,1),
        legend.position="top"
  )
dev.off()

#做检验：
a<-data.frame()
for(i in unique(enh_all_exp_TF_2$tissue)){
  eRNAs_no_p<-(t.test(enh_all_exp_TF_2[enh_all_exp_TF_2$tissue==i&enh_all_exp_TF_2$type=="eRNAs-Enh","TFBS_num"],enh_all_exp_TF_2[enh_all_exp_TF_2$tissue==i&enh_all_exp_TF_2$type=="no_eRNA","TFBS_num"],alternative = "greater")$p.value)
  eRNAs_random_p<-(t.test(enh_all_exp_TF_2[enh_all_exp_TF_2$tissue==i&enh_all_exp_TF_2$type=="eRNAs-Enh","TFBS_num"],enh_all_exp_TF_2[enh_all_exp_TF_2$tissue==i&enh_all_exp_TF_2$type=="random","TFBS_num"],alternative = "greater")$p.value)
  no_random_p<-(t.test(enh_all_exp_TF_2[enh_all_exp_TF_2$tissue==i&enh_all_exp_TF_2$type=="no_eRNA","TFBS_num"],enh_all_exp_TF_2[enh_all_exp_TF_2$tissue==i&enh_all_exp_TF_2$type=="random","TFBS_num"],alternative = "greater")$p.value)
  print(paste(i,"  eRNAs:no_eRNAs:",eRNAs_no_p,"  eRNAs:random:",eRNAs_random_p,"    no-eRNAs:random:",no_random_p,sep=" "))
  a_row<-data.frame(tissue=i,EnheRNA_EnhnoeRNA=eRNAs_no_p,EnheRNA_random=eRNAs_random_p,EnhnoeRNA_random=no_random_p)
  a<-rbind(a,a_row)
}




#做相关性分析，可惜结果太差，不用看了
enh_all_exp_TF$logTPM<-log(enh_all_exp_TF$eRNA_exp_TPM+1,10)
ggplot(enh_all_exp_TF,aes(x=logTPM,y=TFBS_num))+geom_smooth(aes(color=tissue),method = lm,level=0.95)

for(i in unique(enh_all_exp_TF$tissue)){
  cor.value<-cor.test(enh_all_exp_TF[enh_all_exp_TF$tissue==i,"logTPM"],enh_all_exp_TF[enh_all_exp_TF$tissue==i,"TFBS_num"],method = "spearman")$estimate
  cor.p.value<-cor.test(enh_all_exp_TF[enh_all_exp_TF$tissue==i,"logTPM"],enh_all_exp_TF[enh_all_exp_TF$tissue==i,"TFBS_num"],method = "spearman")$p.value
  print(paste(i,cor.value,cor.p.value))
}


# 
# #######################################
# enh_bid_table<-data.frame(TFBS_num=data.frame(t(table(enh_bid_exp_TF$TFBS_num)))$Var2,TFBS_Freq=data.frame(t(table(enh_bid_exp_TF$TFBS_num)))$Freq,type="bid")
# enh_unbid_table<-data.frame(TFBS_num=data.frame(t(table(enh_unbid_exp_TF$TFBS_num)))$Var2,TFBS_Freq=data.frame(t(table(enh_unbid_exp_TF$TFBS_num)))$Freq,type="unbid")
# enh_no_eRNA_table<-data.frame(TFBS_num=data.frame(t(table(enh_no_eRNA_exp_TF$TFBS_num)))$Var2,TFBS_Freq=data.frame(t(table(enh_no_eRNA_exp_TF$TFBS_num)))$Freq,type="no_eRNA")
# 
# enh_all_table<-rbind(enh_bid_table,enh_unbid_table,enh_no_eRNA_table)
# #dcast:
# enh_all_table_dcast<-dcast(enh_all_table,TFBS_num~type,value.var = "TFBS_Freq")
# enh_all_table_dcast[is.na(enh_all_table_dcast)]<-0
# enh_all_table_dcast<-enh_all_table_dcast[order(as.numeric(levels(enh_all_table_dcast$TFBS_num))),]
# #melt:
# enh_all_table_melt<-melt(enh_all_table_dcast,id.vars = "TFBS_num",variable.name = "type",value.name = "TFBS_Freq")
# #bid/unbid/no_eRNA转录因子结合位点数目统计直方图
# ggplot(enh_all_table_melt,aes(x=TFBS_num,y=TFBS_Freq))+geom_bar(stat="identity",position = "fill",aes(fill=type))

#######################################






