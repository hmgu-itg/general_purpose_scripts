#!/usr/bin/env bash

indir=$1
n=$2
t=$3
pheno=$4

indir=${indir/%/}

total=$(ls $indir/*.vcf.gz| wc -l)

if (( $n > $total ));then
    echo "ERROR: there are $total < $n VCF files in $indir" | ts
    exit 1
fi

fname=$(ls $indir/*.vcf.gz| head -n $n| tail -n 1)

process_chunk.sh -i $fname -t $t -p $pheno

