#!/usr/bin/env bash
#
# using bcftools and qctool2 to collect missingness stats from one VCF file
#
###################

set -eo pipefail

indir=$1
phenofile=$2
n=$SLURM_ARRAY_TASK_ID

indir=${indir%/}
total=$(ls $indir/*.vcf.gz| wc -l)

echo "INPUT DIR $indir" | ts
echo "PHENOTYPE FILE $phenofile" | ts
echo "TOTAL FILES $total" | ts
echo "CURRENT FILENO $n" | ts

if (( $n > $total ));then
    echo "ERROR: there are $total < $n VCF files in $indir" | ts
    exit 1
fi

fname=$(ls $indir/*.vcf.gz| sort | head -n $n | tail -n 1)
outname=${fname/%.vcf.gz/.qctool.out}
echo "CURRENT FILE $fname" | ts
echo "OUTPUT FILE $outname" | ts
echo "--------------------------------------------------------------"
echo ""

pheno_name=$(head -n 1 $phenofile| tr ' ' '\t' |cut -f 2)
echo "INFO: phenotype name: $pheno_name"  | ts

bcftools norm -m- "$fname" | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Ov | qctool2 -g - -filetype vcf -differential "$pheno_name" -osnp "$outname" -s "$phenofile"
echo "INFO: done" | ts
echo "--------------------------------------------------------------"
echo ""

if [[ $? -eq 0 ]];then
    echo "INFO: removing $fname" | ts
else
    echo "INFO: something went wrong; keeping $fname" | ts
fi

echo "--------------------------------------------------------------"
echo ""

exit 0
