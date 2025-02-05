#!/bin/Rscript --vanilla

library("data.table")
library("susieR")
library("getopt")

spec<-matrix(c(
  'input','i',1,"character","Input file",
  'exclude','x',1,"character","List of variant IDs to exclude",
  'iterations','t',1,"integer","Max number of iterations (default: 100)",
  'maxnz','m',1,"integer","Max number of non-zero effects (default: 10)",
  'help','h',0,"logical","Help message",
  'lambda','l',0,"logical","Estimate lambda",
  'noprior','n',0,"logical","Skip prior variance check",
  'coverage','c',1,"double","Required coverage (default: 0.95)",
  'purity','p',1,"double","Purity threshold (default: 0.5)"
), byrow=TRUE, ncol=5)
opt<-getopt(spec)

if ( !is.null(opt$help) | is.null(opt$input)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

## IDs to exclude
exclude_list<-character()
if ( ! is.null(opt$exclude) ){
    exclude_list<-tryCatch(read.table(opt$exclude,header=FALSE)$V1, error=function(e) character())
}
cat(sprintf("Variant IDs to exclude: %d\n",length(exclude_list)))

if ( is.null(opt$purity   ) ) { opt$purity   = 0.5    }
if ( is.null(opt$coverage   ) ) { opt$coverage   = 0.95    }
if ( is.null(opt$iterations   ) ) { opt$iterations   = 100    }
if ( is.null(opt$maxnz   ) ) { opt$maxnz   = 10    }
if ( is.null(opt$lambda ) ) { opt$lambda = FALSE }
if ( is.null(opt$noprior ) ) { opt$noprior = FALSE }

## input file
infile<-opt$input
cat(sprintf("Input file: %s\n",infile))

## number of iterations
n_iter<-opt$iterations
cat(sprintf("Max iterations: %d\n",n_iter))

## max L
max_nz<-opt$maxnz
cat(sprintf("Max number of non-zero effects: %d\n",max_nz))

## coverage
C<-opt$coverage
cat(sprintf("Coverage: %.2f\n",C))
cat(sprintf("Check prior variance: %s\n",!opt$noprior))
cat(sprintf("Estimate lambda: %s\n",opt$lambda))

## purity
P<-opt$purity
cat(sprintf("Purity threshold: %.2f\n",P))

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
if (length(z[z==TRUE])>0){
    id_fixed<-df$ID[!z]
    beta_fixed<-df$beta[!z]
    se_fixed<-df$se[!z]
    M<-M[!z,!z]
}else
{
    id_fixed<-df$ID
    beta_fixed<-df$beta
    se_fixed<-df$se
}

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

## excluding variants specified with -x
if (length(exclude_list)!=0){
cat(sprintf("Excluding: %d variants\n",length(exclude_list)))
idx_exclude<-match(exclude_list,id_fixed)
idx_exclude<-idx_exclude[!is.na(idx_exclude)]
if (length(idx_exclude)!=0){
id_fixed<-id_fixed[-idx_exclude]
beta_fixed<-beta_fixed[-idx_exclude]
se_fixed<-se_fixed[-idx_exclude]
M<-M[-idx_exclude,-idx_exclude]
}
}

cat(sprintf("Analyzing: %d variants\n",nrow(M)))

if (opt$lambda){
    lambda<-estimate_s_rss(z=beta_fixed/se_fixed,R=M,n=df$N[1])
    cat(sprintf("Estimated lambda: %.2E\n",lambda))
}

res<-susie_rss(bhat=beta_fixed,shat=se_fixed,n=df$N[1],R=M,max_iter=n_iter,check_prior=!opt$noprior,L=max_nz)

if (res$converged == TRUE){
    cat(paste(c("Convergence:",res$converged),collapse=" "))
    cat(sprintf(" after %d iterations\n",res$niter))
}else
{
    stop("ERROR: failed to converge")
}

## print(res)
cs<-susie_get_cs(res,Xcorr=M,n_purity=nrow(M),coverage=C,min_abs_corr=P)

cat(sprintf("Found %d credible sets\n",length(cs$cs)))
## reporting credible sets and marginal PIPs
## space delimited
## columns are : <"CS"> <CS number> <CS purity> <variant ID> <PIP for this variant in this CS> <marginal PIP for this variant>
## print(cs)
for (idx in 1:length(cs$cs_index)){V<-unlist(cs$cs[idx]);write.table(cbind(rep("CS",length(V)),rep(idx,length(V)),rep(cs$purity$min.abs.corr[idx],length(V)),id_fixed[V],res$alpha[cs$cs_index[idx],V],res$pip[V]),sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE);}
