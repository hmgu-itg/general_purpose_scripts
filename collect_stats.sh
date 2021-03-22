#!/usr/bin/env bash
#
# top level script for collecting missingness stats
#
########################

indir=$1
outdir=$2
pheno=$3
chr=$4

if [[ ! -z "$chr" ]];then
    chroms="$chr"
else
    chroms="1-22"
fi

indir=${indir%/}
outdir=${outdir%/}

echo "INPUT DIR $indir"
echo "OUTPUT DIR $outdir"
echo "PHENOTYPE FILE $pheno"
echo "CHROMS $chroms"
echo "--------------------------------------"
echo ""

logdir="$outdir"/logs
mkdir -p "$logdir"
if [[ ! -d "$outdir" ]];then
    echo "ERROR: could not create output directory $outdir; exit"
    exit 1
fi

tempdir="$outdir"/temp
mkdir -p "$tempdir"
if [[ ! -d "$tempdir" ]];then
    echo "ERROR: could not create temporary directory $tempdir; exit"
    exit 1
fi

ID=$(sbatch --job-name=collect_stats --cpus-per-task=1 --mem-per-cpu=1G --time=10:00:00 -p normal_q --array="$chroms" -o "$logdir"/collect_stats_%A_chr_%a.log -e collect_stats_%A_chr_%a.err /compute/Genomics/software/scripts/general_purpose_scripts/collect_stats_chr.sh "$indir" "$outdir" "$pheno" | cut -d ' ' -f 4)

# cleanup

