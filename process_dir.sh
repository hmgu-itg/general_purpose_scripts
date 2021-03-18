#!/usr/bin/env bash

# current chr
c=$SLURM_ARRAY_TASK_ID

indir=$1
t=$2
pheno=$3

indir=${indir%/}

echo "INPUT DIR=$indir"
echo "P THRESHOLD=$t"
echo "PHENOTYPE FILE=$pheno"
echo "CURRENT CHROM=$n"
echo ""

outdir="$indir"/output_chr"$c"
mkdir "$outdir"
if [[ ! -d "$outdir" ]];then
    echo "ERROR: could not create output directory $outdir"
    exit 1
fi

total=0
for f in $(find "$indir" -maxdepth 1 -mindepth 1 -name "*chr$c*.vcf.gz*");do
    b=$(basename $f)
    ln -s $f "$outdir/$b"
    total=$((total+1))
done

if [[ "$total" -eq 0 ]];then
    echo "INFO: no VCFs for chr $c found in $indir; exit"
    exit 0
fi

echo "-------------------------------------"
echo ""

sbatch -p normal_q --array=1-"$total" -o process_chunk_wrapper_chr"$c"_%A_%a.log -e process_chunk_wrapper_chr"$c"_%A_%a.err /compute/Genomics/software/scripts/general_purpose_scripts/process_chunk_wrapper.sh "$outdir" "$t" "$pheno"

