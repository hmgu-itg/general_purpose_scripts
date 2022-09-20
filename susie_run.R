#!/bin/Rscript --vanilla

library("data.table")
library("susieR")

args<-commandArgs(trailingOnly=T)

# input file
infile<-args[1]

df<-fread(infile,na.strings="nan")
M<-data.matrix(df[,-c("ID","beta","se","N")])
rownames(M)<-NULL
colnames(M)<-NULL

cat(sprintf("Input: %d variants\n",nrow(M)))

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

cat(sprintf("Analyzing: %d variants\n",nrow(M)))

res<-susie_rss(bhat=beta_fixed,shat=se_fixed,n=df$N[1],R=M)
cat(paste(c("Convergence:",res$converged,"\n"),collapse=" "))
cs<-susie_get_cs(res,Xcorr=M,n_purity=nrow(M))

cat(sprintf("Found %d credible sets\n",length(cs$cs_index)))
# for (x in cs$cs){print(length(x));print(x);}
# for (x in cs$cs){print(id_fixed[x])}
i<-1
for (x in cs$cs){write.table(cbind(rep(i,length(x)),id_fixed[x],res$alpha[i,x]),sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE);i<-i+1;}
