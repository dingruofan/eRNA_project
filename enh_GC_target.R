
library(ggplot2)
library(Cairo)
rm(list = ls())

setwd("/home/ding/all_cmd/script/enh_statistics")
#文件来自 /home/ding/all_cmd/script/enh_GC_target.sh
enh_GC_target<-read.table("./enh_GC_target_eRNA_exp",header = F,stringsAsFactors = F,sep="\t")
names(enh_GC_target)<-c("enh_chr","enh_start","enh_end","enh_name","unknown","eRNA_strand","eRNA_exp","gene_strand","gene_id","gene_name","gene_exp","isbid","tissue","type","GC_mean")

#去除unkonwn列：
enh_GC_target$unknown<-NULL

#去除随机位点:
enh_GC_target<-enh_GC_target[enh_GC_target$isbid!="random",]

#双向转录eRNA的表达量取均值：
bid_enh<-unique(enh_GC_target[enh_GC_target$isbid=="bid","enh_name"])
tissues<-unique(enh_GC_target$tissue)
enh_GC_target_1<-enh_GC_target
for(j in tissues){
  for(i in bid_enh){
    # print(i)
    enh_GC_target_1[enh_GC_target_1$enh_name==i&enh_GC_target_1$tissue==j,"eRNA_exp"]<-mean(as.numeric(enh_GC_target_1[enh_GC_target_1$enh_name==i&enh_GC_target_1$tissue==j,"eRNA_exp"]),na.rm = T)
  }
  print(j)
}
#将双向转录enh的eRNA_strand设为"."
enh_GC_target_1[enh_GC_target_1$isbid=="bid","eRNA_strand"]<-"."
#删除重复项:
enh_GC_target_2<-unique(enh_GC_target_1)
enh_GC_target_2$GC_mean<-as.numeric(enh_GC_target_2$GC_mean)
enh_GC_target_2$gene_exp<-as.numeric(enh_GC_target_2$gene_exp)
#做每个组织GC_mean和gene_exp的相关性：除了Placenta。

isbid<-unique(enh_GC_target_2$isbid)

for(i in tissues[-8]){
  for(j in isbid){
    a<-cor.test(enh_GC_target_2[enh_GC_target_2$tissue==i&enh_GC_target_2$isbid==j,"GC_mean"],enh_GC_target_2[enh_GC_target_2$tissue==i&enh_GC_target_2$isbid==j,"gene_exp"],method = "pearson")
    print(paste(i,j,sep="      "))
    print(a)
  }
}

enh_GC_target_3<-enh_GC_target_2[enh_GC_target_2$isbid!="no_eRNA",]
enh_GC_target_3$eRNA_exp<-as.numeric(enh_GC_target_3$eRNA_exp)
enh_GC_target_3$GC_mean<-as.numeric(enh_GC_target_3$GC_mean)
#比较eRNA表达量与GC含量相关性:
is2D_1D<-c("bid","unbid")
for(i in tissues){
  for(j in is2D_1D){
    a<-cor.test(enh_GC_target_3[enh_GC_target_3$tissue==i&enh_GC_target_3$isbid==j,"GC_mean"],enh_GC_target_3[enh_GC_target_3$tissue==i&enh_GC_target_3$isbid==j,"eRNA_exp"])
    print(paste(i,j,sep="      "))
    print(a)
  }
}



png_path<-paste(figure_path,"_histone_DHS.png",sep="")
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 6, height = 6, units='in', dpi=600)

dev.off()
