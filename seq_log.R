#从shell读入工作目录
rm(list=ls())

library(ggplot2)
library(reshape)
library(Cairo)
library(grid)
library(seqLogo)
library(RWebLogo)

source("https://bioconductor.org/biocLite.R")
biocLite("seqLogo")
biocLite("RWebLogo")


