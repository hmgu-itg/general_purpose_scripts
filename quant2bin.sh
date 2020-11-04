#!/bin/bash

source "./functions.sh"

args=("$@")

read -r input frac <<<$(parseCommandLine ${args[@]} "if")

echo $input
echo $frac
frac1=$(echo 1.0-$frac | bc)
echo $frac1

Rscript -e 'library("data.table")
df<-read.table("'$input'",sep="\t",header=T)
f<-function(x,a,b){if (x<a) return(0);if (x <= b & x >= a) return(NA);if (x>b) return(1);}
df$b<-lapply(df$zres,function(x) f(x,quantile(df$zres,'$frac'),quantile(df$zres,'$frac1')))
print(df[c("ID","b")])
'
