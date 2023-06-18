#1.work file
rm(list = ls())
getwd()
setwd("/home/axle/JCMM")
#2.input data
library(data.table)
options(stringsAsFactors = F)
library(tibble)
library(dplyr) 
library(tidyr)
library(readr)
library(tidyverse)
#3.data conversion
#3.1 input gene data
#expr <- fread('fpkm.tsv',data.table = F)
df<-data.table::fread(file='GDC-PANCAN.htseq_fpkm-uq.tsv.gz',data.table = F)

meta<-read.table('gencode.v22.annotation.gene.probeMap',sep='',header = T)
#
meta1 <- meta[1:2]
meta1<-as.data.frame(meta1)
df1<-as.data.frame(df)
df1[1:5,1:5]
rownames(df1) <- df1[,1]        
df1$id<-rownames(df1)
head(df1)
#exp_symbol1<-cbind(df1,meta1)
exp_symbol<-merge(df1,meta1,by="id", all.x=TRUE)
#exp_symbol1 = exp_symbol1[,c(11771:11772,1:11770)]
#rm(exp_symbol1)
exp_symbol[1:5,1:5]
colnames(exp_symbol)
exp_symbol1 = exp_symbol[,c(11771,1:11770)]
exp_symbol1[1:5,1:5]
ids = exp_symbol1[!duplicated(exp_symbol1$gene),]
#write.table(ids,file= "ids.txt")
ids[1:5,1:5]
rownames(ids) <- ids[,1]
rm(df)
rm(df1)
rm(exp_symbol)
rm(exp_symbol1)
library(dplyr)
ids1 <- select(ids,-c(gene,id,xena_sample))
ids1[1:5,1:5]
#4.log，FPKM-TPM
dim(ids1)
expr2=ids1
expr2 <- 2^expr2-1 
dim(expr2)
fpkmToTpm <- function(fpkm){exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}
expr2 <- apply(expr2,2,fpkmToTpm)
expr2[1:5,1:5]
#5.
A <- data.frame(expr2)
write.csv(A,"A-tpm.csv")
#6.TPM-log2（TPM+1）
B <- log2((A)+1)
B[1:5,1:5]
B <- data.frame(B)
write.csv(B,"B-log2(tpm+1).csv")
#7.ZNRF2
select_gene_name<-c("ZNRF2")
B_ZNRF2<-B[select_gene_name,]
#write.table(select_matrix_ZNRF2,file= "select_matrix_ZNRF2.txt")
write.csv(B_ZNRF2,file= "B_ZNRF2.csv")