#!/usr/bin/env bash

indir=$1
indir=${indir%/}

echo "INPUT DIR=$indir"
echo ""

output="$indir"/merged.vcf.gz

total=0
# pattern: _c9_ for chrom 9 etc.
for f in $(find "$indir" -maxdepth 1 -mindepth 1 -name "*_c$c_*.vcf.gz");do
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

sbatch --cpus-per-task=1 --mem-per-cpu=10G --time=10:00:00 -p normal_q --array=1-"$total" -o process_chunk_wrapper_chr"$c"_%A_%a.log -e process_chunk_wrapper_chr"$c"_%A_%a.err /compute/Genomics/software/scripts/general_purpose_scripts/process_chunk_wrapper.sh "$outdir" "$t" "$pheno"

