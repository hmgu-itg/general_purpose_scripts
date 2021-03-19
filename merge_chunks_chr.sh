#!/usr/bin/env bash

indir=$1

indir=${indir%/}
output="$indir"/merged.vcf.gz
plout="$indir"/merged

echo "MERGING VCF FILES IN $indir AND CONVERTING TO PLINK FORMAT"
echo "INPUT DIR=$indir"
echo "OUTPUT=$output"
echo "PLINK OUTPUT PREFIX=$plout"
echo "-------------------------------------"
echo ""

list=$(mktemp -t merge_list-XXXXXXX)
if [[ ! -f $list ]];then
    echo "ERROR: could not create list file; exit"
    exit 1
fi

for f in $(find "$indir" -maxdepth 1 -mindepth 1 -name "*.vcf.gz");do
    realpath $f >> list
done

bcftools concat -f list -Oz -o "$output"
tabix "$output"
rm -f list
plink2 --vcf "$output" --make-pgen --out "$plout"
