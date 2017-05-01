#从shell读入工作目录
rm(list=ls())
setwd("/media/ding/000B49000006264C/eRNA_project/fantom/All_tissue_CAGE")
#############
#双向的：BPKM=TPM/transcript_length
#载入TPM_matrix
tpm_matrix<-read.table("./cage_enh/bidir_pairs_expression_tpm.matrix",header = T,stringsAsFactors = F)
#载入每个组织的ctss的library
# library_count<-read.table("./cage_enh_quantify/library.counts.txt",header = F,stringsAsFactors = F)
# library_count<-as.data.frame(t(library_count))
# names(library_count)<-colnames(tpm_matrix)[2:11]
# rownames(library_count)<-"count"
#载入每个双向转录的双向转录本的长度和:
bid_trans_length<-read.table("./cage_enh/bidir.pairs_transc_length.bed",header = F,stringsAsFactors = F)
names(bid_trans_length)<-c("chr","chr_start","chr_end","genome_location","transcript_length")
bid_trans_length<-subset(bid_trans_length,select = c("genome_location","transcript_length"))

#生成BPKM矩阵:
# for(i in colnames(library_count)){
#   #我草，其实直接有count矩阵，直接除以两个转录本长度和就行:
#   tpm_matrix[,i]<-tpm_matrix[,i]*library_count[,i]/10^6
# }
tpm_matrix<-merge(tpm_matrix,bid_trans_length,by="genome_location")
BPKM_maatrix<-tpm_matrix[,2:11]/tpm_matrix[,12]
BPKM_maatrix$genome_location<-tpm_matrix$genome_location
BPKM_maatrix<-BPKM_maatrix[,c(11,1:10)]

#写入文件:
write.table(BPKM_maatrix,"./cage_enh/bidir_pairs_expression_BPKM.matrix",sep="\t",quote=F,row.names = F,col.names = T)

#单向的:
rm(list=ls())
tpm_matrix<-read.table("./cage_enh/TCs_sorted_exp_tpm.matrix",header = T,stringsAsFactors = F)
#载入单转录的转录本长度和:
unbid_trans_length<-read.table("./cage_enh/TCs_sorted_transc_length.bed",header = F,stringsAsFactors = F)
names(unbid_trans_length)<-c("chr","chr_start","chr_end","genome_location","transcript_length")
unbid_trans_length<-subset(unbid_trans_length,select = c("genome_location","transcript_length"))

tpm_matrix<-merge(tpm_matrix,unbid_trans_length,by="genome_location")
BPKM_maatrix<-tpm_matrix[,2:11]/tpm_matrix[,12]
BPKM_maatrix$genome_location<-tpm_matrix$genome_location
BPKM_maatrix<-BPKM_maatrix[,c(11,1:10)]

write.table(BPKM_maatrix,"./cage_enh/TCs_sorted_exp_BPKM.matrix",sep="\t",quote=F,row.names = F,col.names = T)












