#!/bin/bash

source "./functions.sh"

args=("$@")

read -r input frac <<<$(parseCommandLine ${args[@]} "if")

#echo $input
#echo $frac
#frac1=$(echo 1.0-$frac | bc)
#echo $frac1

Rscript -e '#library("data.table")
df<-read.table("'$input'",sep="\t",header=F,stringsAsFactors=F)
#f<-function(x,a,b){if (x<a) return(0);if (x <= b & x >= a) return(NA);if (x>b) return(1);}
f<-function(x,a){if (x<=a) return(0) else return(1);}
#df$b<-lapply(df$zres,function(x) f(x,quantile(df$zres,'$frac'),quantile(df$zres,'$frac1')))
df$b<-as.numeric(lapply(df$zres,function(x) f(x,quantile(df$zres,'$frac'))))
write.table(df[c("V1","b")],quote=F,sep="\t",row.names=F)
'
