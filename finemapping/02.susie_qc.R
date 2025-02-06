#!/bin/Rscript --vanilla

library("data.table")
library("susieR")
library("ggplot2")
library("getopt")

spec<-matrix(c(
  'input','i',1,"character","Input file",
  'output','o',1,"character","Output prefix",
  'help','h',0,"logical","Help message"
), byrow=TRUE, ncol=5)
opt<-getopt(spec)

if ( !is.null(opt$help) | is.null(opt$input) | is.null(opt$output)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

## input file
infile<-opt$input
cat(sprintf("Input file: %s\n",infile))

## output file
outfile<-opt$output
cat(sprintf("Output prefix: %s\n",outfile))

df<-fread(infile,na.strings="nan")
M<-data.matrix(df[,-c("ID","beta","se","N")])
rownames(M)<-NULL
colnames(M)<-NULL

cat(sprintf("Input: %d variants\n",nrow(M)))

## remove variants that correspond to NA in the correlation matrix
## first remove rows/columns with all entries being NA
z<-rowSums(is.na(M))==ncol(M)
cat(sprintf("Excluding %d variants having all nan:\n",length(z[z==TRUE])))
for (idx in 1:length(z)){
    if (z[idx]){
        cat(sprintf("%s\n",df$ID[idx]))
        }
    }
id_fixed<-df$ID[!z]
beta_fixed<-df$beta[!z]
se_fixed<-df$se[!z]
M<-M[!z,!z]

## remove rows/cols one by one, starting with the one having max number of NAs
z<-rowSums(is.na(M))
idx<-which.max(z)
max_nan<-max(z)
while (max_nan>0){
    cat(sprintf("Excluding %s with %d nan values\n",id_fixed[idx],max_nan))
    id_fixed<-id_fixed[-idx]
    beta_fixed<-beta_fixed[-idx]
    se_fixed<-se_fixed[-idx]
    M<-M[-idx,-idx]
    z<-rowSums(is.na(M))
    max_nan<-max(z)
    idx<-which.max(z)
}

cat(sprintf("Analyzing: %d variants\n",nrow(M)))
zhat<-beta_fixed/se_fixed
lambda<-estimate_s_rss(z=zhat,R=M,n=df$N[1])
cat(sprintf("Estimated lambda: %.2E\n",lambda))
res<-kriging_rss(z=zhat,n=df$N[1],R=M,s=lambda)
tbl<-res$conditional_dist
tbl$ID<-id_fixed
tbl$P2 <- 2*pnorm(abs(tbl$z_std_diff),lower.tail=FALSE)
p2_threshold <- 0.05/nrow(tbl)

## output PNG
pngname<-paste0(outfile,".png")
cat(sprintf("Saving plot in: %s\n",pngname))
png(pngname,width=800,height=800)
print(res$plot)
dev.off()

## full stats table 
tblname<-paste0(outfile,".txt.gz")
cat(sprintf("Saving table in: %s\n",tblname))
#write.table(tbl,file=tblname,sep="\t",row.names=F,quote=F)
fwrite(tbl,file=tblname,sep="\t",row.names=F,quote=F)

## variants to be excluded
df2<-data.frame(ID=subset(tbl,P2<p2_threshold)$ID)
tblname<-paste0(outfile,".excl.txt")
cat(sprintf("Saving %d variants to exclude in: %s\n",nrow(df2),tblname))
write.table(df2,file=tblname,sep="\t",row.names=F,quote=F,col.names=F)
