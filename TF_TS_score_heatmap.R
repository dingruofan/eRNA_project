#统计各个组织的TF.
library(ggplot2)
library(pheatmap)
library(reshape2)
library(Cairo)
rm(list=ls())
setwd("/media/ding/000B49000006264C/eRNA_project/histone")
TF_no_eRNA_count_exp<-read.table("/home/ding/all_cmd/script/enh_statistics/TF_no_eRNA_count_exp",stringsAsFactors = F,header = T)
TF_bid_count_exp<-read.table("/home/ding/all_cmd/script/enh_statistics/TF_bid_count_exp",stringsAsFactors = F,header = T)
TF_unbid_count_exp<-read.table("/home/ding/all_cmd/script/enh_statistics/TF_unbid_count_exp",stringsAsFactors = F,header = T)
TF_count_exp<-merge(merge(TF_no_eRNA_count_exp,TF_bid_count_exp,by=names(TF_no_eRNA_count_exp),all=T),TF_unbid_count_exp,by=names(TF_no_eRNA_count_exp),all=T)
rm(TF_no_eRNA_count_exp,TF_unbid_count_exp,TF_bid_count_exp)

#读入TF_count文件，表示哪些组织里有哪些TF:
TF_no_eRNA_count<-read.table("/home/ding/all_cmd/script/enh_statistics/TF_no_eRNA_count",stringsAsFactors = F,header = T)
TF_bid_count<-read.table("/home/ding/all_cmd/script/enh_statistics/TF_bid_count",stringsAsFactors = F,header = T)
TF_unbid_count<-read.table("/home/ding/all_cmd/script/enh_statistics/TF_unbid_count",stringsAsFactors = F,header = T)
TF_count<-merge(merge(TF_no_eRNA_count,TF_bid_count,by=names(TF_no_eRNA_count),all=T),TF_unbid_count,by=names(TF_no_eRNA_count),all=T)
rm(TF_no_eRNA_count,TF_bid_count,TF_unbid_count)
#需要删除重复的TF:
TF_count_unique<-data.frame(TF_all_tissue=unique(TF_count$TF_all_tissue))
for(i in TF_count_unique$TF_all_tissue){
  for(j in names(TF_count)[c(-1,-length(names(TF_count)))]){
    if(sum(TF_count[TF_count$TF_all_tissue==i,j])>=1){
      TF_count_unique[TF_count_unique$TF_all_tissue==i,j]<-1
    }
    if(sum(TF_count[TF_count$TF_all_tissue==i,j])<1){
      TF_count_unique[TF_count_unique$TF_all_tissue==i,j]<-0
    }
  }
}
#setdiff(TF_count_exp_tmp$TF_all_tissue,TF_count_unique$TF_all_tissue)
#count完全包含exp的TF，这是必须的!

#计算每个转录因子在每个组织中的特异性值
TF_count_exp_tmp<-data.frame(TF_all_tissue=TF_count_exp$TF_all_tissue)
for(j in TF_count_exp$TF_all_tissue ){
  sum_score=0
  sum_exp<-sum(TF_count_exp[TF_count_exp$TF_all_tissue==j,names(TF_count_exp)[2:11]])
  for( i in names(TF_count_exp)[2:11]){
    TSVT<-log(TF_count_exp[TF_count_exp$TF_all_tissue==j,i]/sum_exp,2)
    TF_count_exp_tmp[TF_count_exp_tmp$TF_all_tissue==j,i]<-TSVT
    sum_score<-sum_score+TSVT
  }
  TF_count_exp_tmp[TF_count_exp_tmp$TF_all_tissue==j,"TS_score"]<-sum_score
}
write.table(TF_count_exp_tmp,"/home/ding/all_cmd/script/enh_statistics/TF_TSVT",col.names = T,row.names = F,quote = F,sep="\t")
#TSVT热图：

#ggplot做图:
#总的特异性分值:TSPV，并按照TSPV排序:
tile_melt_TSPV<-TF_count_exp_tmp[,c(1,ncol(TF_count_exp_tmp))]
names(tile_melt_TSPV)[2]<-"TSPV"
tile_melt_TSPV$TF_all_tissue<-factor(tile_melt_TSPV$TF_all_tissue,levels=tile_melt_TSPV[order(tile_melt_TSPV$TSPV),"TF_all_tissue"])
#每个组织的特异性分值:TSVT
tile_melt_TSVT<-melt(TF_count_exp_tmp[,-ncol(TF_count_exp_tmp)],id.vars="TF_all_tissue",variable_name = "tissues",value.name="TSVT")
names(tile_melt_TSVT)[3]<-"TSVT"
tile_melt_TSVT$TF_all_tissue<-factor(tile_melt_TSVT$TF_all_tissue,levels=tile_melt_TSVT[order(tile_melt_TSPV$TSPV,decreasing = T),"TF_all_tissue"])

png_path="/media/ding/000B49000006264C/eRNA_project/figure/TF_TS_score_heatmap_2.png"
CairoPNG(png_path, width = 2, height = 16.7, units='in', dpi=600)

p2<-ggplot(tile_melt_TSPV,aes(x="TSPV",y=reorder(TF_all_tissue,-1*TSPV),fill=TSPV))
p2+geom_tile()+scale_fill_gradient2(low="black",high="white",guide = "colourbar")+  
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.ontop = TRUE)
dev.off()

#这个颜色有点不对:
p1<-ggplot(tile_melt_TSVT,aes(x=tissues,y=TF_all_tissue,fill=TSVT))
p1+geom_tile()+scale_fill_gradient2(midpoint=-7,low="red",mid="green",high="white",guide = "colourbar")+  theme(panel.background = element_blank(),
                                                                                        panel.grid.major.x = element_blank(),
                                                                                        panel.grid.minor.x = element_blank(),
                                                                                        panel.grid.minor.y = element_blank(),
                                                                                        panel.ontop = TRUE)



#做热图：                                                                                                                                                                      panel.ontop = TRUE)
heatmap_TSVT<-TF_count_exp_tmp[order(TF_count_exp_tmp$TS_score),1:12]
heatmap_TSVT_count<-merge(heatmap_TSVT,TF_count_unique,by="TF_all_tissue",suffixes=c(".x",""))
heatmap_TSVT_count<-heatmap_TSVT_count[order(TF_count_exp_tmp$TS_score),c(13:22)]

rownames(heatmap_TSVT)<-heatmap_TSVT$TF_all_tissue
heatmap_TSVT$TF_all_tissue<-NULL
heatmap_TSVT$TS_score<-NULL

pheatmap(heatmap_TSVT,scale = "none",cluster_rows = FALSE,color = colorRampPalette(c("white","green","#FF5151"))(1280),display_numbers = matrix(ifelse(heatmap_TSVT_count >0, "*", ""), nrow(heatmap_TSVT)) )

png_path="/media/ding/000B49000006264C/eRNA_project/figure/TF_TS_score_heatmap.png"
CairoPNG(png_path, width = 16.7, height = 2.2, units='in', dpi=600)

pheatmap(t(heatmap_TSVT),scale = "none",cluster_cols = FALSE,color = colorRampPalette(c("white", "green","#FF5151"))(1280),display_numbers = t(matrix(ifelse(heatmap_TSVT_count >0, "*", ""), nrow(heatmap_TSVT))),fontsize_col = 8)

dev.off()



#查看每个组织最特异的转录因子:
#选择32个特异性的转录因子:TSPV<-40.20171
TF_count_exp_tmp_32<-TF_count_exp_tmp[TF_count_exp_tmp$TS_score<(-40.20171),]


#获得每个组织中TS_score最小的那个转录因子作为这个组织最特异的转录因子
tissues<-names(TF_count_exp_tmp_32)[c(-1,-12)]
a<-data.frame()
for(i in tissues){
  a_1<-TF_count_exp_tmp_32[order( TF_count_exp_tmp_32[,c(i)],decreasing = T),c("TF_all_tissue",i,"TS_score")][1,c(1,3)]
  a_1$tissue<-i
  a_1$TSVT<-TF_count_exp_tmp_32[order( TF_count_exp_tmp_32[,c(i)],decreasing = T),c("TF_all_tissue",i,"TS_score")][1,c(2)]
  a<-rbind(a,a_1)
}

write.table(a,"/home/ding/all_cmd/script/enh_statistics/TS_TF_each_tissue",row.names = F,col.names = T,quote = F,sep="\t")


#获得这32个转录因子在哪个组织中的特异性最高，也就是将32个转录因子按照组织进行分类:
b<-data.frame()
TFs_32<-unique(TF_count_exp_tmp_32$TF_all_tissue)

b1<-melt(TF_count_exp_tmp_32[,-12],id.vars = "TF_all_tissue",variable_name = "tissue")
names(b1)<-c("TF_all_tissue","tissue","TSVT")

for(i in TFs_32){
  b2<-b1[b1$TF_all_tissue==i,]
  b3<-b2[order(b2$TSVT,decreasing=T),][1,]
  b<-rbind(b,b3)
}

b<-(b[order(b$tissue),])

write.table(b,"/home/ding/all_cmd/script/enh_statistics/TS_32TF_each_9tissue",row.names = F,col.names = T,quote = F,sep="\t")






