#!/usr/bin/env bash

n=$SLURM_ARRAY_TASK_ID

indir=$1
t=$2
pheno=$3

indir=${indir%/}
total=$(ls $indir/*.vcf.gz| wc -l)
outdir="$indir"/output

echo "INPUT DIR=$indir"
echo "OUTPUT DIR=$outdir"
echo "P THRESHOLD=$t"
echo "PHENOTYPE FILE=$pheno"
echo "TOTAL FILES=$total"
echo "CURRENT FILENO=$n"

if (( $n > $total ));then
    echo "ERROR: there are $total < $n VCF files in $indir"
    exit 1
fi

fname=$(ls $indir/*.vcf.gz| sort | head -n $n | tail -n 1)

echo "FNAME=$fname"
echo "-------------------------------------"
echo ""

singularity exec -B /compute/Genomics /compute/Genomics/containers/worker_3.1 /compute/Genomics/software/scripts/general_purpose_scripts/process_chunk.sh -i "$fname" -t "$t" -p "$pheno" -o "$outdir"

