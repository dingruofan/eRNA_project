rm(list=ls())
# library(ggplot2)
library(stringr)
library(Cairo)
library(xlsx)
# library(openxlsx)
options(java.parameters = "-Xmx8000m")
library(XLConnect)


#整理enh_all_count_spe.bed:
setwd("/home/ding/all_cmd/script/enh_statistics")
enh_all_count_spe.bed<-read.table("./enh_all_count_spe.bed",header = F,stringsAsFactors = F)
#读入特异性分值文件：
enh_spe<-read.table("./enh_spe",header = T,stringsAsFactors = F)
enh_spe<-enh_spe[,c(1,2,3,4,6)]
names(enh_spe)<-c("chr","chr_start","chr_end","enhancer","Shannon Entropy score")
names(enh_all_count_spe.bed)<-c("chr","chr_start","chr_end","enhancer","has eRNA","tissue specificity","tissue")
#merge合并
enh_all_count_spe.bed<-merge(enh_all_count_spe.bed,enh_spe,by=c("chr","chr_start","chr_end","enhancer"))
enh_all_count_spe.bed<-enh_all_count_spe.bed[,c(1,2,3,4,5,8,6,7)]

enh_all_count_spe.bed$`has eRNA`<-factor(enh_all_count_spe.bed$`has eRNA`,levels = c("bid","unbid","no_eRNA"),labels = c("2D-eRNA","1D-eRNA","no eRNA"))
enh_all_count_spe.bed$`tissue specificity`<-factor(enh_all_count_spe.bed$`tissue specificity`,levels=c("spe","other","uni"),labels = c("Specific","Other","Ubiquitous"))

write.xlsx2(enh_all_count_spe.bed,"/media/ding/000B49000006264C/eRNA_project/figure/table/enhancer_info.xlsx",sheetName = "Enhancer_location",col.names = T,row.names = F,append=F)

#输出为导入数据库的文本格式:
#注意是两个相同数据不同格式的表:
write.table(enh_all_count_spe.bed,"/media/ding/000B49000006264C/eRNA_project/figure/table/database/enhancer_info.tsv",sep = "\t",row.names = F,col.names = T,append = F,quote=F)
#enh_all_count_spe_melt的格式不需要改变，直接拷贝即可。
file.copy("./enh_all_count_spe_melt","/media/ding/000B49000006264C/eRNA_project/figure/table/database/enhancer_info_melt.tsv",overwrite=T)


############################################################
#处理各个组织的./enh_find/enh_all_TF.bed
rm(list=ls())
library(stringr)
# library(xlsx)
library(openxlsx)
options(java.parameters = "-Xmx8000m")
library(XLConnect)


tissues<-str_split("Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney"," ")[[1]]
if(file.exists("/media/ding/000B49000006264C/eRNA_project/figure/table/enhancer_TF_all_tissues.xlsx")){
  file.remove("/media/ding/000B49000006264C/eRNA_project/figure/table/enhancer_TF_all_tissues.xlsx")
}

wb<-createWorkbook()
for(i in tissues){
  setwd(paste("/media/ding/000B49000006264C/eRNA_project/histone/",i,"/enh_find",sep=""))
  enh_all_TF.bed<-read.table("./enh_all_TF.bed",stringsAsFactors = F,header = F,sep="\t")
  names(enh_all_TF.bed)<-c("Enh_chr","Enh_chr_start","Enh_chr_end","enhancer","hsa eRNA","eRNA","TF_chr","TF_chr_start","TF_chr_end","TF","TF_exp")
  
  addWorksheet(wb, sheetName = i, gridLines = T)
  writeData(wb, sheet = i, x = enh_all_TF.bed,
                 colNames = TRUE, rowNames = F)
  #输出为导入数据库的文本格式:
  #因为是追加写入，所以先写入标头:
  enh_all_TF.bed$tissue<-i
  if(!file.exists("/media/ding/000B49000006264C/eRNA_project/figure/table/database/enhancer_TF_all_tissues.tsv")){
    write.table(t(data.frame(names(enh_all_TF.bed))),"/media/ding/000B49000006264C/eRNA_project/figure/table/database/enhancer_TF_all_tissues.tsv",col.names = F,row.names = F,sep = "\t",append = F,quote=F)
  }
  write.table(enh_all_TF.bed,"/media/ding/000B49000006264C/eRNA_project/figure/table/database/enhancer_TF_all_tissues.tsv",sep = "\t",row.names = F,col.names = F,append = T,quote=F)
  
  print(i)
}
saveWorkbook(wb, "/media/ding/000B49000006264C/eRNA_project/figure/table/enhancer_TF_all_tissues.xlsx", overwrite = TRUE)

# #根据上面生成的 enhancer_TF_all_tissues.tsv与enhancer进行合并，添加enh_spe信息:
# rm(list=ls())
# enh_all<-read.table("/home/ding/all_cmd/script/enh_statistics/enh_all_count_spe",stringsAsFactors = F,header = T)
# enh_all<-enh_all[,c("enhancer","type","is_spe")]
# enh_TF<-read.table("/media/ding/000B49000006264C/eRNA_project/figure/table/database/enhancer_TF_all_tissues.tsv",header = T,stringsAsFactors = F,sep="\t")
# enh_TF_merge<-merge(enh_TF,enh_all,by="enhancer")
# #读入random:
# enh_random<-read.table("/home/ding/all_cmd/script/enh_statistics/enh_random_exp_TF_all",header = T,stringsAsFactors = F)
# 
# enh_TF_merge_table<-data.frame(table(enh_TF_merge$enhancer,enh_TF_merge$tissue))
# 
# spe_table<-data.frame(table(enh_all$enhancer))

###########################################
#获得每个组织eRNA的信息：
#来源是enh_eRNA_exp.sh
#2D-eRNA来自每个组织的./enh_find/enh_bid_exp2.bed 
#1D-eRNA来自每个组织的：./enh_find_1/enh_unbid_exp.bed
rm(list=ls())
library(openxlsx)
library(stringr)

#判断xlsx是否存在，存在就删除：
if(file.exists("/media/ding/000B49000006264C/eRNA_project/figure/table/eRNA_info.xlsx")){
  file.remove("/media/ding/000B49000006264C/eRNA_project/figure/table/eRNA_info.xlsx")
}
if(file.exists("/media/ding/000B49000006264C/eRNA_project/figure/table/database/eRNA_info.tsv")){
  file.remove("/media/ding/000B49000006264C/eRNA_project/figure/table/database/eRNA_info.tsv")
}


#写入xlsx:
wb<-createWorkbook()
tissues<-str_split("Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney"," ")[[1]]
ii<-0
for(i in tissues){
  setwd(paste("/media/ding/000B49000006264C/eRNA_project/histone/",i,sep=""))
  # 读入文件：
  bid_exp<-read.table("./enh_find/enh_bid_exp2.bed",header = F,stringsAsFactors = F)
  unbid_exp<-read.table("./enh_find_1/enh_unbid_exp.bed",header = F,stringsAsFactors = F)
  # 补充列名:
  names(unbid_exp)<-c("enh_chr","enh_start","enh_end","enh_name","strand","eRNA_location","dist_to_enh","exp(TPM)")
  names(bid_exp)<-c("enh_chr","enh_start","enh_end","enh_name","strand","eRNA_location","dist_to_enh","exp(TPM)")
  # 补充类型:
  bid_exp$eRNA_type<-"2D-eRNA"
  unbid_exp$eRNA_type<-"1D-eRNA"
  
  #补充2D-eRNA的eRNA的位置信息：
  #/media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/bid_cage_enh_quantify/bidir_pairs_expression_tpm.matrix
  bidir_pairs<-read.table("/media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/cage_enh/bidir_pairs_sorted.bed6_additional",header = F,stringsAsFactors = F)
  names(bidir_pairs)<-c("bid_eRNA_chr","bid_eRNA_start","bid_eRNA_end","eRNA_location","Tag_num","strand")
  bidir_pairs$eRNA_detail<-paste(bidir_pairs$bid_eRNA_chr,":",bidir_pairs$bid_eRNA_start,"-",bidir_pairs$bid_eRNA_end,",",bidir_pairs$strand,sep="")
  #保留双向转录本到最后的eRNA_detail
  bidir_pairs<-bidir_pairs[,c(4,7)]
  cat(nrow(bid_exp))
  bid_exp<-merge(bid_exp,bidir_pairs,by="eRNA_location")
  print(paste(i,nrow(bid_exp)))
  
  #使得单向与双向转录本的格式一致:
  unbid_exp$eRNA_detail<-unbid_exp$eRNA_location
  
  # 合并表:
  eRNA_exp<-rbind(bid_exp,unbid_exp)
  ii<-ii+nrow(eRNA_exp)
  cat(ii)
  cat("\t")
  addWorksheet(wb,sheetName = i,gridLines = T)
  writeData(wb,i,eRNA_exp,colNames = T,rowNames = F)
  
  
  #输出为导入数据库的文本格式:
  eRNA_exp$tisue<-i
  #因为是追加写入，所以先写入标头:
  if(!file.exists("/media/ding/000B49000006264C/eRNA_project/figure/table/database/eRNA_info.tsv")){
    write.table(t(data.frame(names(eRNA_exp))),"/media/ding/000B49000006264C/eRNA_project/figure/table/database/eRNA_info.tsv",col.names = F,row.names = F,sep = "\t",append = F,quote=F)
  }
  write.table(eRNA_exp,"/media/ding/000B49000006264C/eRNA_project/figure/table/database/eRNA_info.tsv",sep = "\t",row.names = F,col.names = F,append = T,quote=F)
  
  print(i)
}
saveWorkbook(wb,"/media/ding/000B49000006264C/eRNA_project/figure/table/eRNA_info.xlsx",overwrite = T)


######################################
#做eRNA在各个组织分布的图：
rm(list=ls())
library(ggplot2)
library(Cairo)
eRNA<-read.table("/media/ding/000B49000006264C/eRNA_project/figure/table/database/eRNA_info.tsv",header =T,stringsAsFactors = F )

eRNA_table<-data.frame(table(eRNA$tisue,eRNA$eRNA_type))
names(eRNA_table)<-c("tissue","type","Freq")

#排序:
eRNA_table$tissue<-reorder(eRNA_table$tissue,eRNA_table$Freq)

#数变做enhancer数目和eRNA数目的相关性：
enh<-read.table("/home/ding/all_cmd/script/enh_statistics/bid_statistics",header = T,stringsAsFactors = F)

#######
enh_bid<-read.table("/home/ding/all_cmd/script/enh_statistics/bid_statistics",header = T,stringsAsFactors = F)
enh_unbid<-read.table("/home/ding/all_cmd/script/enh_statistics/unbid_statistics",header = T,stringsAsFactors = F)
enh_all<-merge(enh_bid,enh_unbid,by=c("tissue","TFBS"))
eRNA_all<-as.data.frame.ts(table(eRNA$tisue,eRNA$eRNA_type))
eRNA_all$tissue<-row.names(eRNA_all)
enh_eRNA_all<-merge(enh_all,eRNA_all,by="tissue")

###########



enh<-enh[,c("tissue","all_enh_count")]

enh_eRNA<-merge(eRNA_table,enh,by="tissue")

#1D-eRNA的数量和增强子数量的相关性:
enh_eRNA_1D<-enh_eRNA[enh_eRNA$type=="1D-eRNA",]
cor.test(enh_eRNA_1D$Freq,enh_eRNA_1D$all_enh_count,method="pearson",alternative = "two.sided")

#2D-eRNA的数量和增强子数量的相关性:
enh_eRNA_2D<-enh_eRNA[enh_eRNA$type=="2D-eRNA",]
cor.test(enh_eRNA_2D$Freq,enh_eRNA_2D$all_enh_count,method="pearson",alternative = "two.sided")

#所有eRNA的数量和增强子数量的相关性:
enh_eRNA_all<-merge(enh_eRNA_1D,enh_eRNA_2D,by=c("tissue","all_enh_count"))
enh_eRNA_all$eRNA_all<-enh_eRNA_all$Freq.x+enh_eRNA_all$Freq.y
cor.test(enh_eRNA_all$eRNA_all,enh_eRNA_all$all_enh_count,method="pearson",alternative = "two.sided")

#1D-eRNA和2D-eRNA的相关性：
cor.test(enh_eRNA_1D$Freq,enh_eRNA_2D$Freq,method="pearson",alternative = "two.sided")

#读入enh_ts信息：
enh_spe<-read.table("/home/ding/all_cmd/script/enh_statistics/enh_all_count_spe_melt",header = T,stringsAsFactors = F)
enh_spe_table<-as.data.frame.ts(table(enh_spe$tissue,enh_spe$is_spe))
enh_spe_table$tissue<-rownames(enh_spe_table)

#获得enh特异性数据：
enh_eRNA_all_ts<-merge(enh_eRNA_all,enh_spe_table,by="tissue")
#2D-eRNA  other:
cor.test(enh_eRNA_all_ts$Freq.y,enh_eRNA_all_ts$other,method="pearson",alternative = "two.sided")
#2D-eRNA  spe:
cor.test(enh_eRNA_all_ts$Freq.y,enh_eRNA_all_ts$spe,method="pearson",alternative = "two.sided")
#2D-eRNA  uni:
cor.test(enh_eRNA_all_ts$Freq.y,enh_eRNA_all_ts$uni,method="pearson",alternative = "two.sided")
#1D-eRNA  other:
cor.test(enh_eRNA_all_ts$Freq.x,enh_eRNA_all_ts$other,method="pearson",alternative = "two.sided")
#1D-eRNA  spe:
cor.test(enh_eRNA_all_ts$Freq.x,enh_eRNA_all_ts$spe,method="pearson",alternative = "two.sided")
#1D-eRNA  uni:
cor.test(enh_eRNA_all_ts$Freq.x,enh_eRNA_all_ts$uni,method="pearson",alternative = "two.sided")
#all-eRNA  other:
cor.test(enh_eRNA_all_ts$eRNA_all,enh_eRNA_all_ts$other,method="pearson",alternative = "two.sided")
#all-eRNA  spe:
cor.test(enh_eRNA_all_ts$eRNA_all,enh_eRNA_all_ts$spe,method="pearson",alternative = "two.sided")
#all-eRNA  uni:
cor.test(enh_eRNA_all_ts$eRNA_all,enh_eRNA_all_ts$uni,method="pearson",alternative = "two.sided")


# 具有相关性 #画图
png_path="/media/ding/000B49000006264C/eRNA_project/figure/eRNA.png"
CairoPNG(png_path, width = 5.7, height = 6, units='in', dpi=600)

ggplot(eRNA_table,aes(tissue,Freq))+geom_bar(stat="identity",aes(fill=type))+theme_bw()+
  ylab("The number of eRNAs")+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_text(size=rel(1.3)),
  axis.text.y = element_text(size=rel(1.3)),
  axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1,size=rel(1.3)),
  legend.text=element_text(size=rel(1.3)),
  legend.title=element_blank(),
  legend.background=element_blank(),
  legend.key = element_blank(),
  legend.margin=margin(0,0,0,0,"mm"),
  legend.box.spacing = unit(1,"mm"),
  #legend.position=c(1,1),legend.justification=c(1,1),
  legend.position="top"
)

dev.off()

######################################


rm(list=ls())
library(stringr)
#获得每个组织eRNA的数据，作为igv的输入的转录本:
tissues<-str_split("Adrenal Brain Breast Heart Liver Lung Ovary Placenta SkeletalMuscle Kidney"," ")[[1]]
for(i in tissues){
  setwd(paste("/media/ding/000B49000006264C/eRNA_project/histone/",i,sep=""))
  # 读入文件：
  bid_exp<-read.table("./enh_find/enh_bid_exp2.bed",header = F,stringsAsFactors = F)
  unbid_exp<-read.table("./enh_find_1/enh_unbid_exp.bed",header = F,stringsAsFactors = F)
  # 补充列名:
  names(unbid_exp)<-c("enh_chr","enh_start","enh_end","enh_name","strand","eRNA_location","dist_to_enh","exp(TPM)")
  names(bid_exp)<-c("enh_chr","enh_start","enh_end","enh_name","strand","eRNA_location","dist_to_enh","exp(TPM)")
  # 补充类型:
  bid_exp$eRNA_type<-"2D-eRNA"
  unbid_exp$eRNA_type<-"1D-eRNA"
  
  #补充2D-eRNA的eRNA的位置信息：
  #/media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/bid_cage_enh_quantify/bidir_pairs_expression_tpm.matrix
  bidir_pairs<-read.table("/media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE/cage_enh/bidir_pairs_sorted.bed6",header = F,stringsAsFactors = F)
  names(bidir_pairs)<-c("eRNA_chr","eRNA_start","eRNA_end","eRNA_location","Tag_num","strand")
  bidir_pairs<-bidir_pairs[,c(4,1,2,3)]
  cat(nrow(bid_exp))
  bid_exp<-merge(bid_exp,bidir_pairs,by="eRNA_location")
  print(paste(i,nrow(bid_exp)))
  
  #使得单向与双向转录本的格式一致:
  unbid_exp$eRNA_chr<-str_split_fixed(unbid_exp$eRNA_location,"[:,-]",4)[,1]
  unbid_exp$eRNA_start<-str_split_fixed(unbid_exp$eRNA_location,"[:,-]",4)[,2]
  unbid_exp$eRNA_end<-str_split_fixed(unbid_exp$eRNA_location,"[:,-]",4)[,3]

  
  # 合并表:
  eRNA_exp<-rbind(bid_exp,unbid_exp)
  
  #获得要输出的列:bed4格式:
  eRNA_exp$score<-"."
  eRNA_exp<-eRNA_exp[,c("eRNA_chr","eRNA_start","eRNA_end","score","score","strand")]
  
  #输出为导入数据库的文本格式:
  write.table(eRNA_exp,paste("/media/ding/000B49000006264C/eRNA_project/figure/table/for_igv/",i,"_eRNA.bed",sep=""),sep = "\t",row.names = F,col.names = F,append = F,quote=F)
  
  print(i)
}








###########################################
#将已经有eRNA注释的eRNA及注释导出为xlsx:
rm(list=ls())
library(openxlsx)
library(stringr)

#原始数据来源:./eRNA_annotate  eRNA_annotate.sh生成。
setwd("/home/ding/all_cmd/script/enh_statistics")
eRNA_ann<-read.table("./eRNA_annotate",stringsAsFactors = F)
names(eRNA_ann)<-c("eRNA_chr","chr_start","chr_end","eRNA_name","unknow","eRNA_strand","gene_chr","gene_start","gene_end","gene_name","gene_type","gene_strand","dist")
eRNA_ann[str_detect(eRNA_ann$eRNA_name,"^bid_eRNA"),"eRNA_name"]<-"2D-eRNA"
eRNA_ann[str_detect(eRNA_ann$eRNA_name,"^unbid_eRNA"),"eRNA_name"]<-"1D-eRNA"
eRNA_ann<-eRNA_ann[,c(-5,-13)]

if(file.exists("/media/ding/000B49000006264C/eRNA_project/figure/table/eRNA_annotated.xlsx")){
  file.remove("/media/ding/000B49000006264C/eRNA_project/figure/table/eRNA_annotated.xlsx")
}

write.xlsx(eRNA_ann,"/media/ding/000B49000006264C/eRNA_project/figure/table/eRNA_annotated.xlsx")


# 2256 /home/ding/all_cmd/script/enh_statistics/bid_eRNA_cluster_transcript.bed
# 8051 /home/ding/all_cmd/script/enh_statistics/unbid_eRNA_cluster_transcript.bed
#10307
# 



########################################################
#获得每个组织增强子的靶基因：
#数据来自enh_target_gene_exp.sh，enh_exp_targetgene_exp
rm(list=ls())
library(stringr)
setwd("/home/ding/all_cmd/script/enh_statistics")
enh_target<-read.table("./enh_exp_targetgene_exp",header = F,stringsAsFactors = F)
names(enh_target)<-c("enh_chr","enh_start","enh_end","enh_name","unknown","eRNA_strand","eRNA_exp(TPM)","gene_strand","gene_id","gene_name","gene_exp(RPKM)","eRNA_type","tissue")
enh_target<-enh_target[,c(-5:-7)]
enh_target<-enh_target[which(str_detect(enh_target$enh_name,"enh_")),]
enh_target<-unique(enh_target)

#载入每个组织每个gene的特异性分值:
# pagenbase<-read.table("/home/ding/all_cmd/script/enh_statistics/pagenbase_all_gene",header = F,stringsAsFactors = F)
# names(pagenbase)<-c("gene_name","TS_score","tissue")

#输出为导入数据库的文本格式:

write.table(enh_target,"/media/ding/000B49000006264C/eRNA_project/figure/table/database/enh_target.tsv",sep = "\t",row.names = F,col.names = T,append = F,quote=F)



#############################
#处理位于CpG岛上的增强子:
rm(list=ls())
library(openxlsx)
CpG_enh<-read.table("/home/ding/all_cmd/script/enh_statistics/CPG_enh_dist0",stringsAsFactors = F,header = F)
names(CpG_enh)<-c("enh_chr","enh_start","enh_end","enh_name","has_eRNA","specificity","tissue","CpG_chr","CpG_start","CpG_end","dist")
#修改内容：
CpG_enh[CpG_enh$has_eRNA=="bid","has_eRNA"]<-"2D-eRNA"
CpG_enh[CpG_enh$has_eRNA=="unbid","has_eRNA"]<-"1D-eRNA"
CpG_enh[CpG_enh$has_eRNA=="no_eRNA","has_eRNA"]<-"no eRNA"
CpG_enh[CpG_enh$specificity=="spe","specificity"]<-"tissue specific"
CpG_enh[CpG_enh$specificity=="other","specificity"]<-"other"
CpG_enh[CpG_enh$specificity=="uni","specificity"]<-"ubiquitous"

CpG_enh<-CpG_enh[,-11]

#写入xlsx表:
if(file.exists("/media/ding/000B49000006264C/eRNA_project/figure/table/CpGi_overlap_enh.xlsx")){
  file.remove("/media/ding/000B49000006264C/eRNA_project/figure/table/CpGi_overlap_enh.xlsx")
}

write.xlsx(CpG_enh,"/media/ding/000B49000006264C/eRNA_project/figure/table/CpGi_overlap_enh.xlsx")

