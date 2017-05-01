#从shell读入工作目录
rm(list=ls())
library(ggplot2)
library(Cairo)
setwd("/home/ding/all_cmd/script/")

eRNA_target<-read.table("./enh_statistics/enh_exp_targetgene_exp",header = F,stringsAsFactors = F)
eRNA_target[eRNA_target$V12=="unbid","V6"]<-eRNA_target[eRNA_target$V12=="unbid","V5"]
names(eRNA_target)<-c("enh_chr","enh_chr_start","enh_chr_end","enh_name","xxx","eRNA_strand","eRNA_exp","gene_strand","gene_id","gene_name","gene_exp","enh_type","tissue")
#先去掉有NA的行:
eRNA_target<-eRNA_target[complete.cases(eRNA_target),]
#对表达量进行归一化:log(RPKM+1,10)
# eRNA_target$gene_exp<-log((eRNA_target$gene_exp+1),10)

#去掉表达量是-1的，也就是Placenta的:
eRNA_target<-subset(eRNA_target,gene_exp>-1)


#每个组织4种enh的target gene的表达量boxplot:
#先排序:
eRNA_target_graph<-eRNA_target
eRNA_target_graph$enh_type<-factor(eRNA_target_graph$enh_type,levels=c("bid","unbid","no_eRNA","random"))

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/target_gene_exp.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width =6.5, height = 6.6, units='in', dpi=600)

ggplot(eRNA_target_graph,aes(x=tissue,y=gene_exp))+
  geom_boxplot(width=0.8,aes(fill=enh_type),show.legend = T,outlier.size = 0.1,outlier.stroke = 1)+
  scale_y_log10(limits=c(10^-3,10^4),breaks=10^(-3:4),labels=10^(-3:4))+
  xlab("tissue")+
  ylab("Expression 0f target genes(RPKM)")+
  scale_fill_discrete(breaks=c("bid","unbid","no_eRNA","random"),labels = c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"]),"random"))+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1.3)),
        axis.text.x=element_text(size=rel(1.3),vjust=1,hjust=1,angle=30),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=rel(1.3)),
        legend.title=element_blank(),
        legend.margin=margin(0,0,0,0,"mm"),
        legend.box.spacing = unit(1,"mm"),
        legend.position="top",
        legend.text=element_text(size=rel(1.3))
  )

dev.off()

#提取各种enh到单个frame,用于P检验：
enh_bid<-eRNA_target[eRNA_target$enh_type=="bid",]
enh_unbid<-eRNA_target[eRNA_target$enh_type=="unbid",]
enh_no_eRNA<-eRNA_target[eRNA_target$enh_type=="no_eRNA",]
enh_random<-eRNA_target[eRNA_target$enh_type=="random",]


#进行P检验
a<-data.frame()
for(i in unique(eRNA_target$tissue)){
  bid_unbid_p<-(t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_type=="bid","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_type=="unbid","gene_exp"],alternative = "greater")$p.value)
  bid_no_eRNA_p<-(t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_type=="bid","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_type=="no_eRNA","gene_exp"],alternative = "greater")$p.value)
  bid_random_p<-(t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_type=="bid","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_type=="random","gene_exp"],alternative = "greater")$p.value)
  unbid_no_eRNA_p<-(t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_type=="unbid","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_type=="no_eRNA","gene_exp"],alternative = "greater")$p.value)
  unbid_random_p<-(t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_type=="unbid","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_type=="random","gene_exp"],alternative = "greater")$p.value)
  no_eRNA_random_p<-(t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_type=="no_eRNA","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_type=="random","gene_exp"],alternative = "greater")$p.value)
  a_row<-data.frame(tissue=i,Enh2D_Enh1D=bid_unbid_p,ENh2D_EnhnoeRNA=bid_no_eRNA_p,ENh2D_random=bid_random_p,Enh1D_EnhnoeRNA=unbid_no_eRNA_p,Enh1D_random=unbid_random_p,EnhnoeRNA_random=no_eRNA_random_p)
  a<-rbind(a,a_row)
  print(paste(i,bid_unbid_p,bid_no_eRNA_p,bid_random_p,unbid_no_eRNA_p,unbid_random_p,no_eRNA_random_p,sep="|"))
}

##################################  
#每个组织eRNA-Enh和no-eRNA Enh以及random进行比较：
eRNA_target$enh_subtype<-""
eRNA_target[eRNA_target$enh_type=="bid"|eRNA_target$enh_type=="unbid","enh_subtype"]<-"eRNA_enh"
eRNA_target[eRNA_target$enh_type=="no_eRNA","enh_subtype"]<-"no_eRNA_enh"
eRNA_target[eRNA_target$enh_type=="random","enh_subtype"]<-"random"
#做图:


png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/target_gene_exp_iseRNA.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 6.5, height = 6.6, units='in', dpi=600)

ggplot(eRNA_target,aes(x=tissue,y=gene_exp))+
  geom_boxplot(width=0.7,aes(fill=enh_subtype),notch=T,outlier.size = 0.5)+
  scale_y_log10(limits=c(10^-3,10^4),breaks=10^(-3:4),labels=10^(-3:4))+
  xlab("tissue")+
  ylab("Expression 0f target genes(RPKM)")+
  scale_fill_discrete(breaks=c("eRNA_enh","no_eRNA_enh","random"),labels=c(expression(Enh["eRNA"]),expression(Enh["no-eRNA"]),"random"))+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1.3)),
        axis.text.x=element_text(size=rel(1.3),vjust=1,hjust=1,angle=30),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=rel(1.3)),
        legend.title=element_blank(),
        legend.margin=margin(0,0,0,0,"cm"),
        legend.position="top",
        legend.text=element_text(size=rel(1.3)),
        legend.box.spacing = unit(0.1,"cm")
  )

dev.off()


#进行P检验
a<-data.frame()
for(i in unique(eRNA_target$tissue)){
  eRNA_no_eRNA.p<-t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_subtype=="eRNA_enh","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_subtype=="no_eRNA_enh","gene_exp"],alternative = "greater")$p.value
  eRNA_random.p<-t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_subtype=="eRNA_enh","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_subtype=="random","gene_exp"],alternative = "greater")$p.value
  no_eRNA_random.p<-t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_subtype=="no_eRNA_enh","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_subtype=="random","gene_exp"],alternative = "greater")$p.value
  a_row<-data.frame(tissue=i,EnheRNA_EnhnoeRNA=eRNA_no_eRNA.p,EnheRNA_random=eRNA_random.p,EnhnoeRNA_random=no_eRNA_random.p)
  a<-rbind(a,a_row)
  print(paste(i,EnheRNA_EnhnoeRNA=eRNA_no_eRNA.p,EnheRNA_random=eRNA_random.p,EnhnoeRNA_random=no_eRNA_random.p))
}







##################################  

#每个组织的enh和random表达量比较:
eRNA_target$enh_subtype<-""
eRNA_target[eRNA_target$enh_type=="bid"|eRNA_target$enh_type=="unbid"|eRNA_target$enh_type=="no_eRNA","enh_subtype"]<-"no_random"
eRNA_target[eRNA_target$enh_type=="random","enh_subtype"]<-"random"
#做图:
eRNA_target$enh_subtype<-factor(eRNA_target$enh_subtype,levels=c("no_random","random"),labels=c("enhancer","random"))

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/target_gene_exp_isEnh.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = par('din')[1], height = par('din')[2], units='in', dpi=600)

ggplot(eRNA_target,aes(x=tissue,y=gene_exp))+
  geom_boxplot(width=0.7,aes(fill=enh_subtype),notch=T,outlier.size = 0.5)+
  scale_y_log10(limits=c(10^-3,10^4),breaks=10^(-3:4),labels=10^(-3:4))+
  xlab("tissue")+
  ylab("Expression 0f target genes(RPKM)")+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1.3)),
        axis.text.x=element_text(size=rel(1.3),vjust=1,hjust=1,angle=30),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=rel(1.3)),
        legend.title=element_blank(),
        legend.margin=margin(0,0,0,0,"mm"),
        legend.box.spacing = unit(1,"mm"),
        legend.position="top",
        legend.text=element_text(size=rel(1.3))
  )
dev.off()

#进行P检验
a<-data.frame()
for(i in unique(eRNA_target$tissue)){
  p.value<-t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_subtype=="enhancer","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$enh_subtype=="random","gene_exp"],alternative = "greater")$p.value
  a_row<-data.frame(tissue=i,enhancer_random=p.value)
  a<-rbind(a,a_row)
  print(paste(i,p.value))
}

  
  # 
# 研究eRNA表达量和gene表达量是否相关：
# ggplot(eRNA_target,aes(x=eRNA_exp,y=gene_exp))+geom_smooth(aes(fill=tissue),method = lm)+scale_y_log10()
# 
# eRNA_target_matrix<-eRNA_target[,c("eRNA_exp","gene_exp")]
# kc<-kmeans(eRNA_target_matrix,2)
# fitted(kc)
# 
# 
# #首先研究eRNA和gene链相同的表达量是否正相关:
# exp_same_strand<-eRNA_target[eRNA_target$eRNA_strand==eRNA_target$gene_strand,]
# ggplot(exp_same_strand,aes(x=eRNA_exp,y=gene_exp))+geom_smooth(aes(fill=tissue),method = lm)
# 
# tissues<-unique(exp_same_strand$tissue)
# 
# for(i in tissues){
#   a<-cor(exp_same_strand[exp_same_strand$tissue==i,"eRNA_exp"],exp_same_strand[exp_same_strand$tissue==i,"gene_exp"])
#   print(a)
# }


######按照是否特异进行分类:

#从shell读入工作目录
rm(list=ls())
library(ggplot2)
library(Cairo)
setwd("/home/ding/all_cmd/script/")

eRNA_target<-read.table("./enh_statistics/enh_exp_targetgene_exp_spe",header = F,stringsAsFactors = F)
eRNA_target[eRNA_target$V12=="unbid","V6"]<-eRNA_target[eRNA_target$V12=="unbid","V5"]
names(eRNA_target)<-c("enh_chr","enh_chr_start","enh_chr_end","enh_name","xxx","eRNA_strand","eRNA_exp","gene_strand","gene_id","gene_name","gene_exp","enh_type","tissue","is_spe")
#先去掉有NA的行:
eRNA_target<-eRNA_target[complete.cases(eRNA_target),]
#对表达量进行归一化:log(RPKM+1,10)
# eRNA_target$gene_exp<-log((eRNA_target$gene_exp+1),10)

#去掉表达量是-1的，也就是Placenta的:
eRNA_target<-subset(eRNA_target,gene_exp>-1)


#每个组织4种enh的target gene的表达量boxplot:
#先排序:
eRNA_target_graph<-eRNA_target
eRNA_target_graph$is_spe<-factor(eRNA_target_graph$is_spe,levels=c("uni","other","spe","random"))

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/target_gene_exp_isspe.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width =6.5, height = 6.6, units='in', dpi=600)

ggplot(eRNA_target_graph,aes(x=tissue,y=gene_exp))+
  geom_boxplot(width=0.7,aes(fill=is_spe),show.legend = T,outlier.size = 0.1,outlier.stroke = 1,size=0.7,notch=F)+
  scale_y_log10(limits=c(10^-3,10^4),breaks=10^(-3:4),labels=10^(-3:4))+
  xlab("tissue")+
  ylab("Expression 0f target genes(RPKM)")+
  scale_fill_discrete(breaks=c("uni","other","spe","random"),labels = c(expression(Enh["UE"]),expression(Enh["Oth"]),expression(Enh["TS"]),"random"))+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1.3)),
        axis.text.x=element_text(size=rel(1.3),vjust=1,hjust=1,angle=30),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=rel(1.3)),
        legend.title=element_blank(),
        legend.margin=margin(0,0,0,0,"mm"),
        legend.position="top",
        legend.text=element_text(size=rel(1.2)),
        legend.box.spacing = unit(0.1,"cm")
  )

dev.off()

#提取各种enh到单个frame,用于P检验：
enh_uni<-eRNA_target[eRNA_target$is_spe=="uni",]
enh_other<-eRNA_target[eRNA_target$is_spe=="other",]
enh_spe<-eRNA_target[eRNA_target$is_spe=="spe",]
enh_random<-eRNA_target[eRNA_target$is_spe=="random",]


#进行P检验
a<-data.frame()
for(i in unique(eRNA_target$tissue)){
  uni_other_p<-(t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$is_spe=="uni","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$is_spe=="other","gene_exp"],alternative = "greater")$p.value)
  uni_spe_p<-(t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$is_spe=="uni","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$is_spe=="spe","gene_exp"],alternative = "greater")$p.value)
  uni_random_p<-(t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$is_spe=="uni","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$is_spe=="random","gene_exp"],alternative = "greater")$p.value)
  other_spe_p<-(t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$is_spe=="other","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$is_spe=="spe","gene_exp"],alternative = "greater")$p.value)
  other_random_p<-(t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$is_spe=="other","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$is_spe=="random","gene_exp"],alternative = "greater")$p.value)
  spe_random_p<-(t.test(eRNA_target[eRNA_target$tissue==i&eRNA_target$is_spe=="spe","gene_exp"],eRNA_target[eRNA_target$tissue==i&eRNA_target$is_spe=="random","gene_exp"],alternative = "greater")$p.value)
  a_row<-data.frame(tissue=i,EnhUE_EnhOth=uni_other_p,ENhUE_EnhTS=uni_spe_p,ENhUE_random=uni_random_p,EnhOth_EnhTS=other_spe_p,EnhOth_random=other_random_p,EnhTS_random=spe_random_p)
  a<-rbind(a,a_row)
  print(paste(i,uni_other_p,uni_spe_p,uni_random_p,other_spe_p,other_random_p,spe_random_p,sep="|"))
}
  

