#!/usr/bin/env bash

rows=$1
cols=$2
N=$3
prefix=$4

for i in $(seq -w 1 $N);do
    echo "Generating file $i/$N" 1>&2
    fname=${prefix}_${i}".txt"
    echo -ne "f.eid\t" > $fname
    for j in $(seq -w $(((i-1)*cols+1)) $((i*cols)));do echo $j;done|perl -lne 'BEGIN{@a=();}{push @a,"COL".$_;}END{print join("\t",@a);}' >> $fname
    for j in $(seq -w 1 $rows | shuf);do
	echo -ne "$j\t" >> $fname
	tr -dc 0-9 </dev/urandom | head -c $(echo "6*$cols" | bc -l) | perl -lne 'print join("\t",/(.{1,6})/g);' >> $fname
    done
done
