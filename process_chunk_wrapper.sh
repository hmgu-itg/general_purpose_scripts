#!/usr/bin/env bash

indir=$1
n=$SLURM_ARRAY_TASK_ID
t=$2
pheno=$3

indir=${indir%/}
total=$(ls $indir/*.vcf.gz| wc -l)

if (( $n > $total ));then
    echo "ERROR: there are $total < $n VCF files in $indir"
    exit 1
fi

fname=$(ls $indir/*.vcf.gz| sort | head -n $n | tail -n 1)

echo "INPUT DIR=$indir"
echo "FNAME=$fname"
echo "P THRESHOLD=$t"
echo "PHENOTYPE FILE=$pheno"

singularity exec -B /compute/Genomics /compute/Genomics/containers/worker_3.1 /compute/Genomics/software/scripts/general_purpose_scripts/process_chunk.sh -i $fname -t $t -p $pheno

