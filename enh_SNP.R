#从shell读入工作目录
rm(list=ls())
library(ggplot2)
library(Cairo)
setwd("/home/ding/all_cmd/script/enh_statistics")

enh_bid_SNP_num<-read.table("enh_bid_cluseter_SNP_num",header = F,stringsAsFactors = F)
enh_unbid_SNP_num<-read.table("enh_unbid_cluseter_SNP_num",header = F,stringsAsFactors = F)
enh_no_eRNA_SNP_num<-read.table("enh_no_eRNA_cluseter_SNP_num",header = F,stringsAsFactors = F)
randomBed_SNP_num<-read.table("./randomBed_cluster_SNP_num",header = F,stringsAsFactors = F)

names(enh_bid_SNP_num)<-c("enhancer","SNP_num")
names(enh_unbid_SNP_num)<-c("enhancer","SNP_num")
names(enh_no_eRNA_SNP_num)<-c("enhancer","SNP_num")
names(randomBed_SNP_num)<-c("enhancer","SNP_num")

enh_bid_SNP_num$type<-"bid"
enh_unbid_SNP_num$type<-"unbid"
enh_no_eRNA_SNP_num$type<-"no_eRNA"
randomBed_SNP_num$type<-"randomBed"

enh_all_SNP_num<-rbind(enh_bid_SNP_num,enh_unbid_SNP_num,enh_no_eRNA_SNP_num,randomBed_SNP_num)
enh_all_SNP_num1<-rbind(enh_bid_SNP_num,enh_unbid_SNP_num,enh_no_eRNA_SNP_num)

#有随机位点的SNP分布情况
# ggplot(enh_all_SNP_num, aes(SNP_num))+ geom_density(aes(fill=type),alpha=0.5)
#无随机位点的SNP分布
# ggplot(enh_all_SNP_num1, aes(SNP_num))+ geom_density(aes(fill=type),alpha=0.5)
enh_all_SNP_num$type<-factor(enh_all_SNP_num$type,levels=c("bid","unbid","no_eRNA","randomBed"))


png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/SNP_isbid_num.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 6, height =5.5, units='in', dpi=600)

ggplot(enh_all_SNP_num,aes(x=type,y=SNP_num))+
  geom_boxplot(aes(fill=type),show.legend = T,outlier.size = 0,outlier.stroke = 0)+
  ylim(0,20)+
  # xlab("distance to enhancer center (bp)")+
  ylab("Number of variants")+
  guides(fill=FALSE)+
  scale_x_discrete(breaks=c("bid","unbid","no_eRNA","randomBed"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"]),"random"))+
  theme_bw()+
  theme(axis.text=element_text(size = rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        axis.title.x=element_blank()
  )
dev.off()

#t.test()
bid_snp<-enh_all_SNP_num[enh_all_SNP_num$type=="bid","SNP_num"]
unbid_snp<-enh_all_SNP_num[enh_all_SNP_num$type=="unbid","SNP_num"]
no_eRNA_snp<-enh_all_SNP_num[enh_all_SNP_num$type=="no_eRNA","SNP_num"]
randomBed_snp<-enh_all_SNP_num[enh_all_SNP_num$type=="randomBed","SNP_num"]

all_type_snp<-c("bid_snp","unbid_snp","no_eRNA_snp","randomBed_snp")

a<-data.frame(stringsAsFactors = F)
i<-1
while(i<length(all_type_snp)){
  j<-i+1
  while(j<=length(all_type_snp)){
    print(paste(all_type_snp[i],all_type_snp[j],t.test(get(all_type_snp[i]),get(all_type_snp[j]),alternative = "greater")$p.value,sep="  "))
    a_row<-data.frame(stringsAsFactors = F,class=paste(all_type_snp[i],all_type_snp[j],sep="_"),p.value=t.test(get(all_type_snp[i]),get(all_type_snp[j]),alternative = "greater")$p.value)
    a<-rbind(a,a_row)
    j<-j+1
  }
  i<-i+1
}
a[a$class=="bid_snp_unbid_snp","class"]<-"Enh2D_Enh1D"
a[a$class=="bid_snp_no_eRNA_snp","class"]<-"Enh2D_EnhnoeRNA"
a[a$class=="bid_snp_randomBed_snp","class"]<-"Enh2D_random"
a[a$class=="unbid_snp_no_eRNA_snp","class"]<-"Enh1D_EnhnoeRNA"
a[a$class=="unbid_snp_randomBed_snp","class"]<-"Enh1D_Enhrandom"
a[a$class=="no_eRNA_snp_randomBed_snp","class"]<-"EnhnoeRNA_random"

#eRNA的和no-eRNA的SNP富集程度比较：
# eRNA_snp<-enh_all_SNP_num[enh_all_SNP_num$type=="bid"|enh_all_SNP_num$type=="unbid","SNP_num"]
# t.test(eRNA_snp,no_eRNA_snp,alternative = "greater")$p.value

############################################################
#做isspe的:
rm(list=ls())
library(ggplot2)
setwd("/home/ding/all_cmd/script/enh_statistics")

enh_spe_SNP_num<-read.table("enh_spe_cluseter_SNP_num",header = F,stringsAsFactors = F)
enh_uni_SNP_num<-read.table("enh_uni_cluseter_SNP_num",header = F,stringsAsFactors = F)
enh_other_SNP_num<-read.table("enh_other_cluseter_SNP_num",header = F,stringsAsFactors = F)
randomBed_SNP_num<-read.table("./randomBed_cluster_SNP_num",header = F,stringsAsFactors = F)

names(enh_spe_SNP_num)<-c("enhancer","SNP_num")
names(enh_uni_SNP_num)<-c("enhancer","SNP_num")
names(enh_other_SNP_num)<-c("enhancer","SNP_num")
names(randomBed_SNP_num)<-c("enhancer","SNP_num")

enh_spe_SNP_num$type<-"spe"
enh_uni_SNP_num$type<-"uni"
enh_other_SNP_num$type<-"other"
randomBed_SNP_num$type<-"randomBed"

enh_all_SNP_num<-rbind(enh_spe_SNP_num,enh_uni_SNP_num,enh_other_SNP_num,randomBed_SNP_num)
enh_all_SNP_num1<-rbind(enh_spe_SNP_num,enh_uni_SNP_num,enh_other_SNP_num)

enh_all_SNP_num$type<-factor(enh_all_SNP_num$type,levels=c("uni","other","spe","randomBed"))


png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/SNP_isspe_num.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 6, height = 5.5, units='in', dpi=600)

ggplot(enh_all_SNP_num,aes(x=type,y=SNP_num))+
  geom_boxplot(aes(fill=type),show.legend = T,outlier.size = 0,outlier.stroke = 0)+
  ylim(0,20)+
  # xlab("distance to enhancer center (bp)")+
  ylab("Number of variants")+
  guides(fill=FALSE)+
  scale_x_discrete(limits=c("uni","other","spe","randomBed"),labels=c(expression(Enh["UE"]),expression(Enh["Oth"]),expression(Enh["TS"],"random")))+
  theme_bw()+
  theme(axis.text=element_text(size = rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        axis.title.x=element_blank()
  )

dev.off()

#做wilcox.test()
spe_snp<-enh_all_SNP_num[enh_all_SNP_num$type=="spe","SNP_num"]
uni_snp<-enh_all_SNP_num[enh_all_SNP_num$type=="uni","SNP_num"]
other_snp<-enh_all_SNP_num[enh_all_SNP_num$type=="other","SNP_num"]
randomBed_snp<-enh_all_SNP_num[enh_all_SNP_num$type=="randomBed","SNP_num"]

all_type_snp<-c("uni_snp","other_snp","spe_snp","randomBed_snp")

i<-1
a<-data.frame(stringsAsFactors = F)
while(i<length(all_type_snp)){
  j<-i+1
  while(j<=length(all_type_snp)){
    print(paste(all_type_snp[i],all_type_snp[j],t.test(get(all_type_snp[i]),get(all_type_snp[j]),alternative = "greater")$p.value,sep="  "))
    a_row<-data.frame(stringsAsFactors = F,p.value=t.test(get(all_type_snp[i]),get(all_type_snp[j]),alternative = "greater")$p.value)
    a<-rbind(a,a_row)
    j<-j+1
  }
  i<-i+1
}

#################################################
#################################################
#做SNP进一步分析：
rm(list=ls())
library(ggplot2)
library(reshape2)
library(stringr)
setwd("/home/ding/all_cmd/script/enh_statistics")
#读入数据:
enh_bid_SNP<-read.table("./enh_bid_cluster_SNP.bed",header = F,stringsAsFactors = F)
enh_unbid_SNP<-read.table("./enh_unbid_cluster_SNP.bed",header = F,stringsAsFactors = F)
enh_no_eRNA_SNP<-read.table("./enh_no_eRNA_cluster_SNP.bed",header = F,stringsAsFactors = F)
enh_random_SNP<-read.table("./randomBed_cluster_SNP.bed",header = F,stringsAsFactors = F)

#重命名:func是Functional category
names(enh_bid_SNP)<-c("name","strand","observed","class","func","enh","dist","type","isspe")
names(enh_unbid_SNP)<-c("name","strand","observed","class","func","enh","dist","type","isspe")
names(enh_no_eRNA_SNP)<-c("name","strand","observed","class","func","enh","dist","type","isspe")
names(enh_random_SNP)<-c("name","strand","observed","class","func","enh","dist","isspe")
#补充enh_random_SNP内容:
enh_random_SNP$type<-"random"

#合并到一个文件：(没有random:)
enh_all_SNP<-rbind(enh_bid_SNP,enh_unbid_SNP,enh_no_eRNA_SNP)

#合并到一个文件：(有random:)
enh_all_SNP_all<-rbind(enh_bid_SNP,enh_unbid_SNP,enh_no_eRNA_SNP,enh_random_SNP)

#
#开始做图：




#先做class分布:(根据isspe)
enh_all_SNP_class<-melt(table(enh_all_SNP$class,enh_all_SNP$isspe))
names(enh_all_SNP_class)<-c("class","isspe","Freq")
enh_all_SNP_class$isspe<-factor(enh_all_SNP_class$isspe,levels=c("spe","other","uni"))
#提取大于0的:
enh_all_SNP_class<-subset(enh_all_SNP_class,Freq>0)
#取log10,那个scale_y_log10的不正确:
# enh_all_SNP_class$Freq<-log(enh_all_SNP_class$Freq,10)

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/SNP_isspe_class.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 6.56, height = 7.3, units='in', dpi=600)

ggplot(enh_all_SNP_class,aes(x=reorder(class,Freq,sum),y=Freq))+
  geom_bar(stat="identity",aes(fill=isspe),width=0.8)+
  scale_y_log10()+
  coord_polar(start=0)+
  scale_x_discrete(breaks=c("single","deletion","insertion","in-del","microsatellite","mnp"),labels=c("single nucleotide variation","deletion","insertion","insertion/deletion","microsatellite","multiple nucleotide polymorphism"))+
  scale_fill_discrete(limits=c("spe","other","uni"),labels=c(expression(Enh["TS"]),expression(Enh["Oth"]),expression(Enh["UE"])))+
  theme_linedraw()+
  theme(axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_text(size=rel(1.2),angle = 20),
        axis.title=element_blank(),
        panel.background=element_rect(fill="white",colour = "white"),
        # panel.border=element_blank(),
        legend.title=element_blank(),
        legend.position="top",
        legend.margin=unit(0.1,"mm"),
        legend.text=element_text(size=rel(1.2))
        )
dev.off()

#再做class分布：(根据type)
enh_all_SNP_type<-melt(table(enh_all_SNP$class,enh_all_SNP$type))
names(enh_all_SNP_type)<-c("class","type","Freq")
enh_all_SNP_type$type<-factor(enh_all_SNP_type$type,levels=c("bid","unbid","no_eRNA"))
#提取大于0的:
enh_all_SNP_type<-subset(enh_all_SNP_type,Freq>0)

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/SNP_isbid_class.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width =6.56, height = 7.3, units='in', dpi=600)

ggplot(enh_all_SNP_type,aes(x=reorder(class,Freq,sum),y=Freq))+
  geom_bar(stat="identity",aes(fill=type),width=0.8)+
  scale_y_log10()+
  coord_polar(theta="x",start=0)+
  scale_x_discrete(breaks=c("single","deletion","insertion","in-del","microsatellite","mnp"),labels=c("single nucleotide variation","Deletion","Insertion","Insertion/Deletion","Microsatellite","Multiple Nucleotide Polymorphism"))+
  scale_fill_discrete(limits=c("bid","unbid","no_eRNA"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"])))+
  theme_linedraw()+
  theme(axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_text(size=rel(1.2),angle = 20),
        axis.title=element_blank(),
        # panel.background=element_rect(fill="white",colour = "white"),
        legend.title=element_blank(),
        legend.position="top",
        legend.margin=unit(0.1,"mm"),
        legend.text=element_text(size=rel(1.2))
  )
dev.off()

#做function(MISO注释的图:）


#做function(MISO注释的图:）
# unknown intergenic
# synonymous_variant	coding-synon
# intron_variant	intron
# downstream_gene_variant	near-gene-3
# upstream_gene_variant	near-gene-5
# nc_transcript_variant	ncRNA
# stop_gained	nonsense
# missense_variant	missense
# stop_lost	stop-loss
# frameshift_variant	frameshift
# inframe_indel	cds-indel
# 3_prime_UTR_variant	untranslated-3
# 5_prime_UTR_variant	untranslated-5
# splice_acceptor_variant	splice-3
# splice_donor_variant	splice-5

#isspe的：
enh_all_SNP$func<-str_split_fixed(enh_all_SNP$func,",",n=5)[,1]
enh_all_SNP_func<-data.frame(t(table(enh_all_SNP$func,enh_all_SNP$isspe)))
names(enh_all_SNP_func)<-c("isspe","func","Freq")
enh_all_SNP_func$func<-as.character(enh_all_SNP_func$func)
enh_all_SNP_func<-enh_all_SNP_func[which(enh_all_SNP_func$Freq>1),]
# enh_all_SNP_func$Freq<-log(enh_all_SNP_func$Freq,10)
#排序：
enh_all_SNP_func$isspe<-factor(enh_all_SNP_func$isspe,levels=c("uni","other","spe"))

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/SNP_isspe_function.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 6.7, height = 6.6, units='in', dpi=600)

ggplot(enh_all_SNP_func,aes(x=reorder(func,Freq,max),y=Freq))+
  geom_bar(position="dodge",stat="identity",aes(fill=isspe),width = 0.7)+
  scale_y_log10(breaks=10^(1:4),labels=10^(1:4))+
  coord_flip()+
  ylab("SNP number")+
  xlab("Function")+
  scale_fill_discrete(limits=c("spe","other","uni"),labels=c(expression(Enh["TS"]),expression(Enh["Oth"]),expression(Enh["UE"])))+
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(1.1)),
        legend.title=element_blank(),
        legend.position="top",
        legend.text=element_text(size=rel(1.3)),
        legend.margin=unit(0.1,"mm")
  )
dev.off()

#isbid的:
enh_all_SNP$func<-str_split_fixed(enh_all_SNP$func,",",n=5)[,1]
enh_all_SNP_func<-data.frame(t(table(enh_all_SNP$func,enh_all_SNP$type)))
names(enh_all_SNP_func)<-c("type","func","Freq")
enh_all_SNP_func$func<-as.character(enh_all_SNP_func$func)
enh_all_SNP_func<-enh_all_SNP_func[which(enh_all_SNP_func$Freq>1),]
# enh_all_SNP_func$Freq<-log(enh_all_SNP_func$Freq,10)

enh_all_SNP_func$type<-factor(enh_all_SNP_func$type,levels=c("bid","unbid","no_eRNA"))

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/SNP_isbid_function.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width =6.7, height = 6.6, units='in', dpi=600)

ggplot(enh_all_SNP_func,aes(x=reorder(func,Freq),y=Freq))+
  geom_bar(position="dodge",stat="identity",aes(fill=type),width = 0.7)+
  coord_flip()+
  scale_y_log10(breaks=10^(1:4),labels=10^(1:4))+
  ylab("SNP number")+
  xlab("Function")+
  scale_fill_discrete(limits=c("bid","unbid","no_eRNA"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"])))+
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(1.1)),
        axis.text.y=element_text(size=rel(1.1)),
        legend.title=element_blank(),
        legend.position="top",
        legend.text=element_text(size=rel(1.3)),
        legend.margin=unit(0.1,"mm")
  )
dev.off()



############################
#做距离enh的图:
rm(list=ls())
library(ggplot2)
library(stringr)
library(Cairo)
setwd("/home/ding/all_cmd/script/enh_statistics")
enh_all_SNP_dist<-read.table("./enh_all_SNP_dist")

names(enh_all_SNP_dist)<-c("SNP_class","SNP_function","type","isspe","dist")

#取function的第一个字段:
enh_all_SNP_dist$SNP_function<-str_split_fixed(enh_all_SNP_dist$SNP_function,",",5)[,1]

#做class图:
png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/SNP_histogram_class.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width =6.7, height = 6.6, units='in', dpi=600)

enh_all_SNP_dist$SNP_class<-factor(enh_all_SNP_dist$SNP_class,levels=c("single","insertion","deletion","in-del","microsatellite","mnp"))

enh_all_SNP_dist_1<-subset(enh_all_SNP_dist,SNP_class=="single"|SNP_class=="deletion"|SNP_class=="insertion")
binwidth<-diff(range(enh_all_SNP_dist_1$dist))/1000

ggplot(enh_all_SNP_dist_1,aes(dist))+
  geom_histogram(position = "identity",aes(fill=SNP_class),binwidth = binwidth,alpha=0.4)+
  xlim(-2000,2000)+
  xlab("Distance to enhancer center (bp)")+
  ylab("Number of variants")+
  theme_bw()+
  scale_fill_manual(breaks=c("single","insertion","deletion"),labels=c("SNV","Insertion","Deletion"),values=c("red","blue","green"),name="ding")+
  # scale_x_continuous(breaks=seq(-2000,2000,1000),labels=paste(seq(-2000,2000,1000),"bp",seq=""),limits = c(-2000,2000))+
  theme(
        axis.text=element_text(size=rel(1.3)),
        axis.title = element_text(size=rel(1.3)),
        # axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.margin=margin(0,0,0,0,"mm"),
        legend.position=c(1,1),
        legend.justification=c(1.26,1.2),
        legend.text.align = 0,
        legend.background=element_blank(),
        legend.text=element_text(size=rel(1.15))
  )

dev.off()


class_table<-data.frame(table(enh_all_SNP_dist$SNP_class,enh_all_SNP_dist$dist))
names(class_table)<-c("SNP_class","dist","Freq")

ggplot(class_table,aes(x=dist,y=Freq))+
  geom_point(aes(color=SNP_class))+xlim(-2000,2000)

#做function图:
png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/SNP_histogram_function.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width =6.7, height = 6.6, units='in', dpi=600)

binwidth<-diff(range(enh_all_SNP_dist$dist))/1000
ggplot(enh_all_SNP_dist,aes(dist))+
  geom_histogram(aes(fill=SNP_function),binwidth = binwidth,alpha=0.7)+
  xlim(-2000,2000)+
  xlab("distance to enhancer center (bp)")+
  ylab("SNP number")+
  scale_fill_discrete(palette(gray(seq(0,.9,len = 25))),breaks=c("unknown","near-gene-5","near-gene-3","intron"),labels=c("unknown","near-gene-5","near-gene-3","intron"))+
  # scale_x_continuous(breaks=seq(-2000,2000,1000),labels=paste(seq(-2000,2000,1000),"bp",seq=""),limits = c(-2000,2000))+
  theme_bw()+
  theme(axis.text=element_text(size=rel(1.1)),
        # axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.margin=unit(0,"mm"),
        legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.background=element_blank(),
        legend.text=element_text(size=rel(1.0))
  )

dev.off()


#画图:画每个组织的SNP数量差异：
rm(list=ls())
library(ggplot2)
library(Cairo)
setwd("/home/ding/all_cmd/script/enh_statistics")
#数据来源于enh_SNP.sh
enh_SNP<-read.table("/media/ding/000B49000006264C/eRNA_project/figure/table/database/enh_all_SNP",header = T,stringsAsFactors = F)
randomBed_SNP_num<-read.table("./randomBed_cluster_SNP_num",header = F,stringsAsFactors = F)

#处理随机位点:
names(randomBed_SNP_num)<-c("enh_name","Freq")
randomBed_SNP_num$tissue<-"random"

#获得所有增强子的数据:
enh_all<-read.table("enh_all_count_spe_melt",header = T,stringsAsFactors = F)
enh_all_1<-enh_all[,c("enhancer","tissue"),drop=F]
names(enh_all_1)<-c("enh_name","tissue")

enh_Freq<-data.frame(table(enh_SNP$enh_name))
names(enh_Freq)<-c("enh_name","Freq")

#合并没有SNP的，记为0:
enh_all_Freq<-merge(enh_all_1,enh_Freq,by="enh_name",all.x=T)
enh_all_Freq[is.na(enh_all_Freq$Freq),"Freq"]<-0

#合并random的:
enh_all_Freq_random<-rbind(enh_all_Freq,randomBed_SNP_num)

#标注组织和random不同的颜色:
enh_all_Freq_random[enh_all_Freq_random$tissue=="random","colours"]<-"grey"
enh_all_Freq_random[enh_all_Freq_random$tissue!="random","colours"]<-"blue"


#做boxplot:
png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/SNP_boxplot.png"
CairoPNG(png_path, width = 5.5, height =5.16 , units='in', dpi=600)

ggplot(enh_all_Freq_random,aes(x=reorder(tissue,Freq,median),y=Freq))+
  geom_boxplot(aes(fill=colours),show.legend = T,outlier.size = 0,outlier.stroke = 0)+
  ylim(0,15)+
  xlab("tissue")+
  ylab("SNP number")+
  scale_fill_manual(values = c("#FF5151","grey"))+
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
  )+
  guides(fill=F)

dev.off()


###后期补充做图，对于先前无影响:
rm(list=ls())
library(ggplot2)
setwd("/home/ding/all_cmd/script/enh_statistics")
bid_SNP<-read.table("enh_bid_cluster_SNP.bed",header = F,stringsAsFactors = F,sep = "\t")
unbid_SNP<-read.table("enh_unbid_cluster_SNP.bed",header = F,stringsAsFactors = F,sep = "\t")
no_eRNA_SNP<-read.table("enh_no_eRNA_cluster_SNP.bed",header = F,stringsAsFactors = F,sep = "\t")
random_SNP<-read.table("randomBed_cluster_SNP.bed",header = F,stringsAsFactors = F,sep = "\t")
#给random添加一列:
random_SNP$v9<-"random"

#重命名:
names(bid_SNP)<-c("SNP_name","strand","obs","SNP_class","function","enh","Freq","isbid","isspe")
names(unbid_SNP)<-c("SNP_name","strand","obs","SNP_class","function","enh","Freq","isbid","isspe")
names(no_eRNA_SNP)<-c("SNP_name","strand","obs","SNP_class","function","enh","Freq","isbid","isspe")
names(random_SNP)<-c("SNP_name","strand","obs","SNP_class","function","enh","Freq","isbid","isspe")


enh_all_SNP_num<-rbind(bid_SNP,unbid_SNP,no_eRNA_SNP,random_SNP)

#有随机位点的SNP分布情况
enh_all_SNP_num$type<-factor(enh_all_SNP_num$type,levels=c("bid","unbid","no_eRNA","randomBed"))

table(enh_all_SNP_num$enh,enh_all_SNP_num$SNP_class,enh_all_SNP_num$isbid)

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/SNP_isbid_num_1.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 6, height =5.5, units='in', dpi=600)

ggplot(enh_all_SNP_num,aes(x=isbid,y=Freq))+
  geom_boxplot(aes(fill=SNP_class),show.legend = T,outlier.size = 0,outlier.stroke = 0)+
  ylim(0,20)+
  # xlab("distance to enhancer center (bp)")+
  ylab("The number of SNPs")+
  guides(fill=FALSE)+
  scale_x_discrete(breaks=c("bid","unbid","no_eRNA","randomBed"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"]),"random"))+
  theme_bw()+
  theme(axis.text=element_text(size = rel(1.3)),
        axis.title=element_text(size=rel(1.3)),
        axis.title.x=element_blank()
  )

dev.off()



#补充做图：
#从shell读入工作目录
rm(list=ls())
library(ggplot2)
library(reshape)
library(Cairo)
setwd("/home/ding/all_cmd/script/enh_statistics")

enh_bid_SNP_num<-read.table("enh_bid_cluseter_SNP_num_detail",header = F,stringsAsFactors = F)
enh_unbid_SNP_num<-read.table("enh_unbid_cluseter_SNP_num_detail",header = F,stringsAsFactors = F)
enh_no_eRNA_SNP_num<-read.table("enh_no_eRNA_cluseter_SNP_num_detail",header = F,stringsAsFactors = F)
randomBed_SNP_num<-read.table("./randomBed_cluster_SNP_num_detail",header = F,stringsAsFactors = F)

names(enh_bid_SNP_num)<-c("enhancer","SNV","Insertion","Deletion")
names(enh_unbid_SNP_num)<-c("enhancer","SNV","Insertion","Deletion")
names(enh_no_eRNA_SNP_num)<-c("enhancer","SNV","Insertion","Deletion")
names(randomBed_SNP_num)<-c("enhancer","SNV","Insertion","Deletion")

enh_bid_SNP_num$type<-"bid"
enh_unbid_SNP_num$type<-"unbid"
enh_no_eRNA_SNP_num$type<-"no_eRNA"
randomBed_SNP_num$type<-"randomBed"

#有random:
enh_all_SNP_num<-rbind(enh_bid_SNP_num,enh_unbid_SNP_num,enh_no_eRNA_SNP_num,randomBed_SNP_num)
#无random:
enh_all_SNP_num1<-rbind(enh_bid_SNP_num,enh_unbid_SNP_num,enh_no_eRNA_SNP_num)
#进行变长:
enh_all_SNP_num_melt<-melt(enh_all_SNP_num,id = c("enhancer","type"))
names(enh_all_SNP_num_melt)[3]<-"variants"

#排序:
enh_all_SNP_num_melt$type<-factor(enh_all_SNP_num_melt$type,levels=c("bid","unbid","no_eRNA","randomBed"))

png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/SNP_isbid_num_detail.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 6, height =5.5, units='in', dpi=600)

ggplot(enh_all_SNP_num_melt,aes(x=type,y=value))+
  geom_boxplot(aes(fill=variants),show.legend = T,outlier.size = 0,outlier.stroke = 0)+
  scale_y_log10(breaks=10^(1:2),labels=10^(1:2))+
  # xlab("distance to enhancer center (bp)")+
  ylab("The number of variants")+
  #guides(fill=FALSE)+
  scale_x_discrete(breaks=c("bid","unbid","no_eRNA","randomBed"),labels=c(expression(Enh["2D-eRNA"]),expression(Enh["1D-eRNA"]),expression(Enh["no-eRNA"]),"random"))+
  theme_bw()+
  theme(
       axis.text.x=element_text(size=rel(1.1),vjust=1,hjust=1,angle=0),
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


############################################################
#补充做图:
#做isspe的:
rm(list=ls())
library(ggplot2)
library(reshape)
setwd("/home/ding/all_cmd/script/enh_statistics")

enh_spe_SNP_num<-read.table("enh_spe_cluseter_SNP_num_detail",header = F,stringsAsFactors = F)
enh_uni_SNP_num<-read.table("enh_uni_cluseter_SNP_num_detail",header = F,stringsAsFactors = F)
enh_other_SNP_num<-read.table("enh_other_cluseter_SNP_num_detail",header = F,stringsAsFactors = F)
randomBed_SNP_num<-read.table("./randomBed_cluster_SNP_num_detail",header = F,stringsAsFactors = F)

names(enh_spe_SNP_num)<-c("enhancer","SNV","Insertion","Deletion")
names(enh_uni_SNP_num)<-c("enhancer","SNV","Insertion","Deletion")
names(enh_other_SNP_num)<-c("enhancer","SNV","Insertion","Deletion")
names(randomBed_SNP_num)<-c("enhancer","SNV","Insertion","Deletion")

enh_spe_SNP_num$type<-"spe"
enh_uni_SNP_num$type<-"uni"
enh_other_SNP_num$type<-"other"
randomBed_SNP_num$type<-"randomBed"

enh_all_SNP_num<-rbind(enh_spe_SNP_num,enh_uni_SNP_num,enh_other_SNP_num,randomBed_SNP_num)
enh_all_SNP_num1<-rbind(enh_spe_SNP_num,enh_uni_SNP_num,enh_other_SNP_num)

enh_all_SNP_num1_melt<-melt(enh_all_SNP_num1,id = c("enhancer","type"))
names(enh_all_SNP_num1_melt)[3]<-"variants"

sum(enh_all_SNP_num1_melt[enh_all_SNP_num1_melt$variants=="SNV","value"])/sum(enh_all_SNP_num1_melt[,"value"])

sum(enh_all_SNP_num1_melt[enh_all_SNP_num1_melt$variants=="Insertion","value"])/sum(enh_all_SNP_num1_melt[,"value"])

sum(enh_all_SNP_num1_melt[enh_all_SNP_num1_melt$variants=="Deletion","value"])/sum(enh_all_SNP_num1_melt[,"value"])


#进行变长:
enh_all_SNP_num_melt<-melt(enh_all_SNP_num,id = c("enhancer","type"))
names(enh_all_SNP_num_melt)[3]<-"variants"

enh_all_SNP_num_melt$type<-factor(enh_all_SNP_num_melt$type,levels=c("uni","other","spe","randomBed"))


png_path<-"/media/ding/000B49000006264C/eRNA_project/figure/SNP_isspe_num_detail.png"
# 使用Cairo生成高清图片：par('din')获得当前g
CairoPNG(png_path, width = 6, height = 5.5, units='in', dpi=600)

ggplot(enh_all_SNP_num_melt,aes(x=type,y=value))+
  geom_boxplot(aes(fill=variants),show.legend = T,outlier.size = 0,outlier.stroke = 0)+
  scale_y_log10(breaks=10^(1:2),labels=10^(1:2))+
  # xlab("distance to enhancer center (bp)")+
  ylab("The number of variants")+
  #guides(fill=FALSE)+
  scale_x_discrete(breaks=c("uni","other","spe","randomBed"),labels=c(expression(Enh["TS"]),expression(Enh["UE"]),expression(Enh["Oth"]),"random"))+
  theme_bw()+
  theme(
    axis.text.x=element_text(size=rel(1.1),vjust=1,hjust=1,angle=0),
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



