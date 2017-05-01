rm(list=ls())
library(rGREAT)
library(reshape2)
library(ggplot2)
setwd("/home/ding/all_cmd/script/enh_statistics")
#读入所有enhancer信息：
enh_all<-read.table("./enh_all_count",header = T,stringsAsFactors = F)

#进行宽变长:
enh_all_melt<-melt(enh_all,id.vars = c("chr","chr_start","chr_end","enhancer","tissue_statistics","type"),)
enh_all_melt<-enh_all_melt[enh_all_melt$value!=0,]
names(enh_all_melt)[c(7,8)]<-c("tissue","num")

#对每个组织进行GO分析:
tissues<-as.character(unique(enh_all_melt$tissue))
for(i in tissues){
  assign(paste(i,"_enh",sep=""),enh_all_melt[enh_all_melt$tissue==i,])
  print(i)
}
#从rGREAT官网建立任务:
# for(j in tissues){
#   assign(paste(j,"_job",sep=""),submitGreatJob(get(paste(j,"_enh",sep=""))[,1:3],request_interval = 5))
#   print(j)
# }


Adrenal_job<-submitGreatJob(Adrenal_enh[,1:3])
Brain_job<-submitGreatJob(Brain_enh[,1:3])
Breast_job<-submitGreatJob(Breast_enh[,1:3])
Heart_job<-submitGreatJob(Heart_enh[,1:3])
Liver_job<-submitGreatJob(Liver_enh[,1:3])
Lung_job<-submitGreatJob(Lung_enh[,1:3])
Ovary_job<-submitGreatJob(Ovary_enh[,1:3])
Placenta_job<-submitGreatJob(Placenta_enh[,1:3])
SkeletalMuscle_job<-submitGreatJob(SkeletalMuscle_enh[,1:3])
Kidney_job<-submitGreatJob(Kidney_enh[,1:3])

#提取enrichment前5个做图：

deal_job<-function(tissue_job){
  tissue_tb<-getEnrichmentTables(tissue_job,ontology = "GO Molecular Function")[[1]][,c(1,2,3,8)]
  tissue_tb$q.value=p.adjust(tissue_tb$Binom_Raw_PValue,method = "fdr")
  tissue_tb$enrich<-(-1)*log(tissue_tb$q.value,10)
  tissue_tb<-tissue_tb[order(tissue_tb$enrich,decreasing = T),][1:5,]
  tissue_tb
}

for(j in tissues){
  tissue_tb<-deal_job(get(paste(j,"_job",sep="")))
  png(paste("./GO/",j,"_GO",sep=""),width = 1000,height=500)
  print(ggplot(tissue_tb,aes(y=enrich,x=name))+geom_bar(stat="identity",width = 0.5)+coord_flip()+theme(axis.title.x = element_text(size=20,colour = "black",face = "bold"),axis.text.y = element_text(size=20,colour = "black"),title=element_text(size=25,colour = "black"))+labs(title = j))
  dev.off()
  print(paste("完成了",j))
}




#做各个组织靶基因的GO分析:

rm(list=ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(RDAVIDWebService)
library(openxlsx)

setwd("/home/ding/all_cmd/script/enh_statistics/")
#如果使用TS enhancer的靶基因:
target_all<-read.table("./target_gene/TS_gene/all_TS_target",header = T,stringsAsFactors = F)
#如果使用所有增强子的靶基因：
# target_all<-read.table("./target_gene/all_gene/all_target",header = T,stringsAsFactors = F)

names(target_all)<-c("SYMBOL","tissue")

entrezID<-unique(bitr(target_all$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))

#转entrezID:
target_all_ID<-merge(target_all,entrezID,by="SYMBOL")

#如果要写入的文件存在，则不写入。
if(file.exists("/media/ding/000B49000006264C/eRNA_project/figure/table/GO_BP_TS_gene.xlsx")){
  file.remove("/media/ding/000B49000006264C/eRNA_project/figure/table/GO_BP_TS_gene.xlsx")
}
if(file.exists("/media/ding/000B49000006264C/eRNA_project/figure/table/GO_MF_TS_gene.xlsx")){
  file.remove("/media/ding/000B49000006264C/eRNA_project/figure/table/GO_MF_TS_gene.xlsx")
}
if(file.exists("/media/ding/000B49000006264C/eRNA_project/figure/table/GO_CC_TS_gene.xlsx")){
  file.remove("/media/ding/000B49000006264C/eRNA_project/figure/table/GO_CC_TS_gene.xlsx")
}



#并写入xlsx表:
wb_BP<-createWorkbook()
wb_MF<-createWorkbook()
wb_CC<-createWorkbook()
for(i in unique(target_all_ID$tissue)){
  
  GO_ont<-c("MF","BP","CC")
  for(j in GO_ont){
    ego <- enrichGO(gene          = target_all_ID[target_all_ID$tissue==i,"ENTREZID"],
                    # universe      = names(geneList),
                    OrgDb         = org.Hs.eg.db,
                    ont           = j,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    # if(nrow(summary(ego))>0){
      assign(paste(i,"_",j,sep=""),summary(ego))
      print(paste(i,"中",j,"有",nrow(summary(ego)),"条",sep=""))
    # }
  }
  #将十个组织BP的写入文件：
  if(nrow(get(paste(i,"_BP",sep="")))>0){
    addWorksheet(wb_BP, sheetName = i, gridLines = T)
    writeData(wb_BP, sheet = i, x = get(paste(i,"_BP",sep="")),colNames = TRUE, rowNames = F)
  }
  if(exists(paste(i,"_BP",sep=""))){
    rm(list=(paste(i,"_BP",sep="")))
  }
  
  
  #将10个组织的MF写入文件：
  if(nrow(get(paste(i,"_MF",sep="")))>0){
    addWorksheet(wb_MF, sheetName = i, gridLines = T)
    writeData(wb_MF, sheet = i, x = get(paste(i,"_MF",sep="")),colNames = TRUE, rowNames = F)
  }
  if(exists(paste(i,"_MF",sep=""))){
    rm(list=paste(i,"_MF",sep=""))
  }
  
  #将10个组织的CC写入文件：
  if(nrow(get(paste(i,"_CC",sep="")))>0){
    addWorksheet(wb_CC, sheetName = i, gridLines = T)
    writeData(wb_CC, sheet = i, x = get(paste(i,"_CC",sep="")),colNames = TRUE, rowNames = F)
  }
  if(exists(paste(i,"_CC",sep=""))){
    rm(list=paste(i,"_CC",sep=""))
  }
  
}
saveWorkbook(wb_BP, "/media/ding/000B49000006264C/eRNA_project/figure/table/GO_BP_TS_gene.xlsx", overwrite = TRUE)
saveWorkbook(wb_MF, "/media/ding/000B49000006264C/eRNA_project/figure/table/GO_MF_TS_gene.xlsx", overwrite = TRUE)
saveWorkbook(wb_CC, "/media/ding/000B49000006264C/eRNA_project/figure/table/GO_CC_TS_gene.xlsx", overwrite = TRUE)

 

#做KEGG富集分析:
#并写入xlsx表:
if(file.exists("/media/ding/000B49000006264C/eRNA_project/figure/table/KEGG_all_gene.xlsx")){
  file.remove("/media/ding/000B49000006264C/eRNA_project/figure/table/KEGG_all_gene.xlsx")
}

wb_kegg<-createWorkbook()
for(i in unique(target_all_ID$tissue)){
  
  kk <- enrichKEGG(gene         = target_all_ID[target_all_ID$tissue==i,"ENTREZID"],
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)
  # head(summary(kk))
  
  if(nrow(summary(kk))>0){
    addWorksheet(wb_kegg, sheetName = i, gridLines = T)
    writeData(wb_kegg, sheet = i, x = summary(kk),colNames = TRUE, rowNames = F)
  }
  if(exists("kk")){
    rm(list=c("kk"))
  }
  
}
saveWorkbook(wb_kegg, "/media/ding/000B49000006264C/eRNA_project/figure/table/KEGG_all_gene.xlsx", overwrite = TRUE)


####################
#做每个组织靶基因的挑出来的组织相关的top10 GO分析图:
rm(list=ls())
library(openxlsx)
library(ggplot2)
library(Cairo)
setwd("/home/ding/all_cmd")

GO_top10<-read.xlsx("./stastics_human.xlsx",sheet="target_GO",startRow = 1,colNames = T,rowNames = F)
GO_top10$point_shape<-"0"
GO_top10$point_size<-"15"
tissues<-unique(GO_top10$tissue)
#每个组织加个最大p值和最小p值:



for(i in tissues){
  
  # 使用Cairo生成高清图片：par('din')获得当前g
  # png_path<-paste("/media/ding/000B49000006264C/eRNA_project/figure/GO_BP/","colourbar.png",sep="")
  
  i=tissues[8]
  png_path<-paste("/media/ding/000B49000006264C/eRNA_project/figure/GO_BP/",i,"_G0_BP.png",sep="")

  CairoPNG(png_path, width =par('din')[1], height = par('din')[2], units='in', dpi=600)
  #hight :3.4
  
  ggplot(GO_top10[GO_top10$tissue==i,],aes(x=0,y=reorder(Description,sort(pvalue,decreasing=T))))+
    geom_point(aes(colour=pvalue,size=point_size,shape=point_shape))+
    scale_shape_manual(values = 15)+
    scale_size_manual(values=8)+
    scale_y_discrete(position = "right") +
    scale_x_continuous(breaks=NULL,label=NULL,limits = c(0,0))+
    # scale_colour_continuous(high = 'lightblue', low = 'darkblue',limits=c(-19,-4),breaks=seq(-19,-4,5),label=c(expression(10^-19),expression(10^-14),expression(10^-9),expression(10^-4)))+
    scale_colour_continuous(high = 'lightblue', low = 'darkblue',limits=10^c(-19,-3),guide="legend")+
    guides(size=F,shape=F,colour=F)+
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=rel(1.5)),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      panel.border = element_blank(),
      panel.spacing=unit(0,"cm"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", 
                                      colour = NA),
      legend.position="top"
    )
  
  # assign(paste(i,"_plotp",sep=""),p)
  dev.off()
  
}


#批量做图
dev.new()
for(i in tissues){
  png_path<-paste("/media/ding/000B49000006264C/eRNA_project/figure/GO_BP/",i,"G0_BP.png",sep="")
  # 使用Cairo生成高清图片：par('din')获得当前g
  CairoPNG(png_path, width =4.21, height = 3.01, units='in', dpi=600)
  
  get(paste(i,"_plotp",sep=""))
  
  dev.off()
}

