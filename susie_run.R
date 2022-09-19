#!/bin/Rscript --vanilla

library("data.table")
library("susieR")

args<-commandArgs(trailingOnly=T)

infile<-args[1]

df<-fread(infile,na.strings="nan")
M<-data.matrix(df[,-c("ID","beta","se","N")])
rownames(M)<-NULL
colnames(M)<-NULL

paste(c("Input:",nrow(M),"variants"),collapse=" ")

# remove variants that correspond to NA in the correlation matrix
# first remove rows/columns with all entries being NA
z<-rowSums(is.na(M))==ncol(M)
id_fixed<-df$ID[!z]
beta_fixed<-df$beta[!z]
se_fixed<-df$se[!z]
M<-M[!z,!z]

# keep rows/columns where there are no NA
z<-rowSums(is.na(M))==0
id_fixed<-id_fixed[z]
beta_fixed<-beta_fixed[z]
se_fixed<-se_fixed[z]
M<-M[z,z]

paste(c("Analyzing:",nrow(M),"variants"),collapse=" ")

res<-susie_rss(bhat=beta_fixed,shat=se_fixed,n=df$N[1],R=M)
paste(c("Convergence:",res$converged),collapse=" ")
cs<-susie_get_cs(res,Xcorr=M)

paste(c("Found",length(cs$cs_index),"credible sets"),collapse=" ")
for (x in cs$cs){print(length(x));print(x);}
