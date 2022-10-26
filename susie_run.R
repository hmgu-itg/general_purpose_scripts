#!/bin/Rscript --vanilla

library("data.table")
library("susieR")
library("getopt")

spec<-matrix(c(
  'input','i',1,"character","Input file",
  'iterations','t',1,"integer","Max number of iterations (default: 100)",
  'help','h',0,"logical","Help message",
  'lambda','l',0,"logical","Estimate lambda"
), byrow=TRUE, ncol=5)
opt<-getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if ( !is.null(opt$help) | is.null(opt$input)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

## set some reasonable defaults for the options that are needed,
## but were not specified.
if ( is.null(opt$iterations   ) ) { opt$iterations   = 100    }
if ( is.null(opt$lambda ) ) { opt$lambda = FALSE }

## args<-commandArgs(trailingOnly=T)

## input file
infile<-opts$input

## number of iterations
n_iter<-opt$iterations
## if (length(args)>1){n_iter<-as.integer(args[2]);}
cat(sprintf("Iterations: %d\n",n_iter))

df<-fread(infile,na.strings="nan")
M<-data.matrix(df[,-c("ID","beta","se","N")])
rownames(M)<-NULL
colnames(M)<-NULL

cat(sprintf("Input: %d variants\n",nrow(M)))

## remove variants that correspond to NA in the correlation matrix
## first remove rows/columns with all entries being NA
z<-rowSums(is.na(M))==ncol(M)
id_fixed<-df$ID[!z]
beta_fixed<-df$beta[!z]
se_fixed<-df$se[!z]
M<-M[!z,!z]

## keep rows/columns where there are no NA
z<-rowSums(is.na(M))==0
id_fixed<-id_fixed[z]
beta_fixed<-beta_fixed[z]
se_fixed<-se_fixed[z]
M<-M[z,z]

cat(sprintf("Analyzing: %d variants\n",nrow(M)))

if (opts$lambda){
    lambda<-estimate_s_rss(z=beta_fixed/se_fixed,R=M,n=df$N[1])
    cat(sprintf("Estimated lambda: %.2E\n",lambda))
}

res<-susie_rss(bhat=beta_fixed,shat=se_fixed,n=df$N[1],R=M,max_iter=n_iter)
if (res$converged == TRUE){
    cat(paste(c("Convergence:",res$converged),collapse=" "))
    cat(sprintf(" after %d iterations\n",res$niter))
}else
{
    stop("ERROR: failed to converge")
}

## print(res)
cs<-susie_get_cs(res,Xcorr=M,n_purity=nrow(M))

cat(sprintf("Found %d credible sets\n",length(cs$cs)))
## reporting credible sets and marginal PIPs
## space delimited
## columns are : <"CS"> <CS number> <CS purity> <variant ID> <PIP for this variant in this CS> <marginal PIP for this variant>
## print(cs)
for (idx in 1:length(cs$cs_index)){V<-unlist(cs$cs[idx]);write.table(cbind(rep("CS",length(V)),rep(idx,length(V)),rep(cs$purity$min.abs.corr[idx],length(V)),id_fixed[V],res$alpha[cs$cs_index[idx],V],res$pip[V]),sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE);}
