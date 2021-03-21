#!/usr/bin/env bash

# current chr
c=$SLURM_ARRAY_TASK_ID

indir=$1
outdir=$2
pheno=$3

indir=${indir%/}
outdir=${outdir%/}

echo "INPUT DIR $indir"
echo "ROOT OUTPUT DIR $outdir"
echo "PHENOTYPE FILE $pheno"
echo "CURRENT CHR $c"

# should exist already
logdir="$outdir"/logs

if [[ ! -d "$logdir" ]];then
    echo "ERROR: log dir $logdir does not exist"
    exit 1
fi

if [[ ! -d "$outdir" ]];then
    echo "ERROR: root output dir $outdir does not exist"
    exit 1
fi

outdir2="$outdir"/chr"$c"
mkdir -p "$outdir2"
if [[ ! -d "$outdir2" ]];then
    echo "ERROR: could not create output directory $outdir2 for chr $c"
    exit 1
fi

echo "OUTPUT DIR $outdir2"
echo "LOG DIR $logdir"
echo "START" $(date)
echo "---------------------------------------------------------------"
echo ""

total=0
# pattern: _c9_ for chrom 9 etc.
for f in $(find "$indir" -maxdepth 1 -mindepth 1 -name "*_c${c}_*.vcf.gz");do
    b=$(basename $f)
    ln -s $f "${outdir2}/$b"
    total=$((total+1))
done

if [[ "$total" -eq 0 ]];then
    echo "INFO: no VCFs for chr $c found in $indir; exit"
    exit 0
fi

sbatch --job-name=merge_chr"$c" --cpus-per-task=1 --mem-per-cpu=10G --time=10:00:00 -p normal_q --array=1-"$total" -o "$logdir"/collect_stats_chr_"$c"_%A_part_%a.log -e "$logdir"/collect_stats_chr_"$c"_%A_part_%a.err --wrap="singularity exec -B /compute/Genomics /compute/Genomics/containers/worker_3.1 /compute/Genomics/software/scripts/general_purpose_scripts/collect_stats_chunk.sh $outdir2 $pheno"

