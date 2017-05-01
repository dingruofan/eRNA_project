#从shell读入工作目录
rm(list=ls())
library(RWebLogo)

setwd("/home/ding/all_cmd/script/enh_statistics/meme_motif")

library(Biostrings)
data(HNF4alpha)
pfm <- consensusMatrix(HNF4alpha)
640*4
D
