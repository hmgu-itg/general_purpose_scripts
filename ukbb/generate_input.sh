#!/usr/bin/env bash

rows=$1
cols=$2
N=$3
prefix=$4

m=$((cols*N))
n=${#m}

for i in $(seq -w 1 $N);do
    echo "Generating file $i/$N" 1>&2
    fname=${prefix}_${i}".txt.gz"
    echo -ne "f.eid\t" | gzip -c > $fname
    for j in $(seq $(((i-1)*cols+1)) $((i*cols)));do echo $j;done|perl -slne 'BEGIN{@a=();}{push @a,sprintf("COL%0".$w."d",$_);}END{print join("\t",@a);}' -- -w=$n | gzip -c >> $fname
    for j in $(seq -w 1 $rows | shuf);do
	echo -ne "$j\t"
	tr -dc 0-9 </dev/urandom | head -c $(echo "6*$cols" | bc -l)|perl -lne 'print join("\t",/(.{1,6})/g);'|gawk -v seed=$RANDOM 'BEGIN{FS=OFS="\t";srand(seed);}{x=rand();n=1+int(NF*rand());if (x>=0 && x<0.05){$n="";} if (x>=0.05 && x<0.1){$n="\""$n"\"";} if(x>=0.1 && x<0.15){$n=rand();} if(x>=0.15 && x<0.2){$n=rand()", "rand();} if(x>=0.2 && x<0.25){$n="\""rand()", "rand()"\"";} if(x>=0.25 && x<0.3){$n=" ";} print $0;}'
    done | gzip - -c >> $fname
done
