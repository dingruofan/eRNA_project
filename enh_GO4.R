rm(list=ls())
library(rGREAT)
library(reshape2)
library(ggplot2)
library(scales)
library(Cairo)
setwd("/home/ding/all_cmd/script/enh_statistics")
#读入所有enhancer信息：
enh_all<-read.table("./enh_all_count_spe_melt",header = T,stringsAsFactors = F)

#进行宽变长:


#对每个组织spe的enh进行GO分析,首先分配变量并赋予data.frame
tissues<-as.character(unique(enh_all$tissue))
for(i in tissues){
  assign(paste(i,"_enh",sep=""),enh_all[enh_all$tissue==i&enh_all$is_spe=="uni",])
  print(i)
}
#从rGREAT官网建立任务:
for(j in tissues){
  assign(paste(j,"_job",sep=""),submitGreatJob(get(paste(j,"_enh",sep=""))[,1:3],request_interval = 5))
  print(j)
}


#从这里开始执行：
#提取enrichment

deal_job<-function(tissue_job,tissue){
  tissue_tb<-getEnrichmentTables(tissue_job,ontology = "GO Molecular Function")[[1]][,c(1,2,3,8)]
  tissue_tb$q.value=p.adjust(tissue_tb$Binom_Raw_PValue,method = "fdr")
  tissue_tb$enrich<-(-1)*log(tissue_tb$q.value,10)
  tissue_tb<-subset(tissue_tb[order(tissue_tb$enrich,decreasing = T),],q.value<0.01)
  tissue_tb$tissue<-tissue
  tissue_tb
}
# print(ggplot(tissue_tb,aes(y=enrich,x=name))+geom_bar(stat="identity",width = 0.5)+coord_flip()+theme(axis.title.x = element_text(size=20,colour = "black",face = "bold"),axis.text.y = element_text(size=20,colour = "black"),title=element_text(size=25,colour = "black"))+labs(title = j))

#获得每个组织的GO item:
for(j in tissues){
  assign(paste(j,"_tissue_tb",sep=""),deal_job(get(paste(j,"_job",sep="")),j)[,])
}

#合并十个组织的GO：
all_tissue_tb<-data.frame()
for(j in tissues){
  all_tissue_tb<-rbind(all_tissue_tb,get(paste(j,"_tissue_tb",sep="")))
}
all_tissue_tb_GO<-data.frame(t(table(all_tissue_tb$name)))[,c(2,3)]
names(all_tissue_tb_GO)<-c("name","Freq")
all_tissue_tb_GO<-all_tissue_tb_GO[order(all_tissue_tb_GO$Freq,decreasing=TRUE),]

#选择其中Freq>=7的：
all_tissue_GO_common<-all_tissue_tb_GO[all_tissue_tb_GO$Freq>=7,"name",drop=FALSE]

#有两条太长了，去掉。。。
all_tissue_GO_common<-all_tissue_GO_common[-(grep("oxidoreductase activity",all_tissue_GO_common$name)),"name",drop=FALSE]


#获得每个组织common的：
for(j in tissues){
  assign(paste(j,"_tissue_tb_common",sep=""),merge(get(paste(j,"_tissue_tb",sep="")),all_tissue_GO_common,by="name"))
}

#合并十个组织的common GO:
all_tissue_tb_common<-data.frame()
for(j in tissues){
  all_tissue_tb_common<-rbind(all_tissue_tb_common,get(paste(j,"_tissue_tb_common",sep="")))
}

#做图：
png_path="/media/ding/000B49000006264C/eRNA_project/figure/GO_uni_enh_common.png"
CairoPNG(png_path, width = par('din')[1], height = par('din')[2], units='in', dpi=600)

ggplot(all_tissue_tb_common,aes(x=reorder(tissue,enrich,decreasing=TRUE),y=reorder(name,Binom_Genome_Fraction)))+
  geom_point(aes(colour=enrich,size=Binom_Genome_Fraction))+
  scale_colour_gradient2("Enrichment",low="white",mid = "blue", high = "red",guide="colorbar")+
  scale_size("Genome Fraction",breaks=seq(0.2,0.8,0.2))+
  xlab("tissue")+
  ylab("GO name")+
  theme(axis.text=element_text(size=rel(1.1)),
        legend.text=element_text(size=rel(1.1)),
        axis.text=element_text(size = rel(1.1)),
        axis.text.x=element_text(angle = 40,vjust=1,hjust=1),
        axis.title=element_blank(),
        legend.text=element_text(size=rel(1.2)),
        legend.background=element_blank(),
        legend.key = element_blank(),
        legend.margin=unit(0,"mm")
  )
dev.off()




# source("https://bioconductor.org/biocLite.R")
# biocLite("GOFunction")


#Adrenal:NADPH-adrenodoxin reductase activity  NADPH-肾上腺皮质氧还蛋白还原酶活性
#Brain: brain-derived neurotrophic factor binding 脑源性神经营养因子结合
#       brain-derived neurotrophic factor-activated receptor activity 脑源性神经营养因子激活受体活性
#Breast fibrinogen binding 
#heart : oltage-gated potassium channel activity involved in ventricular cardiac muscle cell action potential repolarization,电压门控钾通道活动参与心室心肌细胞动作电位复极化
#Kidney  glucocorticoid receptor binding  糖皮质激素受体结合
#        insulin-like growth factor binding  胰岛素样生长因子受体结合
#        insulin-like growth factor receptor binding  胰岛素样生长因子受体结合
#Liver  
# L-leucine transaminase activity   L-亮氨酸转氨酶活性
# L-valine transaminase activity  L-缬氨酸转氨酶活性
# L-isoleucine transaminase activity  L-异亮氨酸转氨酶活性
# L-tyrosine:2-oxoglutarate aminotransferase activity  L-酪氨酸：2-氧戊二酸氨基转移酶活性
# [heparan sulfate]-glucosamine 3-sulfotransferase 3 activity [硫酸乙酰肝素] - 葡萄糖胺3-磺基转移酶3活性
# transaminase activity  转氨酶活性

#Lung:
# fibroblast growth factor-activated receptor activity 成纤维细胞生长因子激活受体活性

#Ovary:heparin binding 肝素结合 ，可能与卵巢早衰有关，打肝素可以保胎?

#Placenta 

#SkeletalMuscle
# structural constituent of muscle 肌肉的结构成分
# tropomyosin binding 原肌球蛋白结合
# actin binding  肌动蛋白结合
# 以及其他各种与肌醇的关系，再次不列出，太多了


