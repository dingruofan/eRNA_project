#统计各个组织的TF.
rm(list=ls())
setwd("/media/ding/000B49000006264C/eRNA_project/histone")
TF_bid_count_exp<-read.table("/home/ding/all_cmd/script/enh_statistics/TF_bid_count_exp",stringsAsFactors = F,header = T)

TF_bid_count_exp_tmp<-data.frame(TF_all_tissue=TF_bid_count_exp$TF_all_tissue)
for(j in TF_bid_count_exp$TF_all_tissue ){
  sum_score=0
  sum_exp<-sum(TF_bid_count_exp[TF_bid_count_exp$TF_all_tissue==j,names(TF_bid_count_exp)[2:11]])
  for( i in names(TF_bid_count_exp)[2:11]){
    TSVT<-log(TF_bid_count_exp[TF_bid_count_exp$TF_all_tissue==j,i]/sum_exp,2)
    sum_score<-sum_score+TSVT
  }
  TF_bid_count_exp_tmp[TF_bid_count_exp_tmp$TF_all_tissue==j,"TS_score"]<-sum_score
}

TF_bid_count_exp<-cbind(TF_bid_count_exp,TF_bid_count_exp_tmp[,"TS_score",drop=FALSE])

write.table(TF_bid_count_exp,"/home/ding/all_cmd/script/enh_statistics/TF_bid_count_exp",col.names = T,row.names = F,sep = "\t",quote = F)

