#!/usr/bin/env bash

indir=$1
threads=$2

indir=${indir%/}
output="$indir"/merged.vcf.gz
plout="$indir"/merged

echo "MERGING VCF FILES AND CONVERTING TO PLINK FORMAT"
echo "INPUT DIR $indir"
echo "USING THREADS $threads"
echo "OUTPUT VCF $output"
echo "PLINK OUTPUT PREFIX $plout"
echo "-------------------------------------"
echo ""

flist=$(mktemp -t merge_list-XXXXXXX)
if [[ ! -f $flist ]];then
    echo "ERROR: could not create merging list file; exit"
    exit 1
fi

for f in $(find "$indir" -maxdepth 1 -mindepth 1 -name "*.vcf.gz");do
    realpath $f >> "$flist"
done

singularity exec -B /compute/Genomics /compute/Genomics/containers/worker_3.1 bcftools concat -a -f "$flist" -Oz -o "$output" --threads "$threads"
singularity exec -B /compute/Genomics /compute/Genomics/containers/worker_3.1 tabix "$output"
rm -f "$flist"
singularity exec -B /compute/Genomics /compute/Genomics/containers/worker_3.1 plink2 --vcf "$output" --make-pgen --out "$plout" --threads "$threads"
