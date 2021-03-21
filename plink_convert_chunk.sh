#!/usr/bin/env bash
#
# converts one chunk VCF to PGEN
#
################################

indir=$1
threads=$2
fileno=$SLURM_ARRAY_TASK_ID

indir=${indir%/}

echo "CONVERTING FILE $fileno TO PLINK FORMAT"
echo "INPUT DIR $indir"
echo "USING THREADS $threads"

# total number of VCF files
total=$(find "$indir" -maxdepth 1 -mindepth 1 -name "*.vcf.gz" ! -name "merged*"| wc -l)

if (( $fileno > $total ));then
    echo "ERROR: there are $total < $fileno VCF files in $indir"
    exit 1
fi

# file to process
fname=$(find "$indir" -maxdepth 1 -mindepth 1 -name "*.vcf.gz" ! -name "merged*"| sort | head -n $fileno| tail -n 1)
outprefix=${fname/%.vcf.gz}

echo "FILE NO $fileno"
echo "FILE NAME $fname"
echo "START" $(date)
echo "-------------------------------------"
echo ""

# TODO: change --vcf-half-call mode if necessary
singularity exec -B /compute/Genomics /compute/Genomics/containers/worker_3.1 plink2 --vcf "$fname" --make-pgen --out "$outprefix" --threads "$threads" --vcf-half-call m

echo "DONE" $(date)
