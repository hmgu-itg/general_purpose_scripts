#!/bin/bash

source "./functions.sh"

args=("$@")

read -r input output <<<$(parseCommandLine ${args[@]} "io")

echo $input
echo $output


#function(x,a,b){if (x<a) return(0);if (x <= b & x >= a) return(NA);if (x>b) return(1);}
#> df$b<-lapply(df$zres,function(x) f(x,quantile(df$zres,0.45),quantile(df$zres,0.55)))
