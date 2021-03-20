#!/usr/bin/env bash

# current chr
c=$SLURM_ARRAY_TASK_ID

indir=$1
t=$2
pheno=$3
threads=$4

indir=${indir%/}

echo "INPUT DIR $indir"
echo "P THRESHOLD $t"
echo "THREADS $threads"
echo "PHENOTYPE FILE $pheno"
echo "CURRENT CHROM $c"
echo ""

# symlinks go here
outdir="$indir"/output_chr"$c"

# processed chunks go here
outdir2="$outdir"/processed

mkdir -p "$outdir"
if [[ ! -d "$outdir" ]];then
    echo "ERROR: could not create output directory $outdir"
    exit 1
fi

total=0
# pattern: _c9_ for chrom 9 etc.
for f in $(find "$indir" -maxdepth 1 -mindepth 1 -name "*_c${c}_*.vcf.gz");do
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

ID=$(sbatch --job-name=process_chr"$c" --cpus-per-task=1 --mem-per-cpu=10G --time=10:00:00 -p normal_q --array=1-"$total" -o process_chunk_wrapper_chr"$c"_%A_%a.log -e process_chunk_wrapper_chr"$c"_%A_%a.err /compute/Genomics/software/scripts/general_purpose_scripts/process_chunk_wrapper.sh "$outdir" "$t" "$pheno" "$outdir2"| cut -d ' ' -f 4)

ID=$(sbatch --job-name=merge_chr"$c" --dependency=afterok:$ID --cpus-per-task=${threads} --mem=200G --time=10:00:00 -p normal_q --array=1 -o merge_chr"$c"_%A_%a.log -e merge_chr"$c"_%A_%a.err /compute/Genomics/software/scripts/general_purpose_scripts/merge_chunks_chr.sh "$outdir2" "$threads"| cut -d ' ' -f 4)

# cleanup


