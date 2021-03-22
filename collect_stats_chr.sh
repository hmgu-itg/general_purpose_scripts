#!/usr/bin/env bash

function usage {
    echo ""
    echo "Usage: $0 -i <input dir>"
    echo "          -o <output dir>"
    echo "          -p <pheno file>"
    echo "        { -r : if analyses should be reusmed; optional; default: false }"
    exit 0
}

resume="no"
OPTIND=1
while getopts "i:o:p:r" optname; do
    case "$optname" in
        "i" ) input="${OPTARG}";;
        "o" ) output="${OPTARG}";;
        "p" ) pheno="${OPTARG}";;
        "r" ) resume="yes";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

input=${input%/}
output=${output%/}

# current chr
c=$SLURM_ARRAY_TASK_ID

echo "INPUT DIR $input"
echo "ROOT OUTPUT DIR $output"
echo "PHENOTYPE FILE $pheno"
echo "CURRENT CHR $c"

# should exist already
logdir="$output"/logs
tempdir="$output"/temp

if [[ ! -d "$logdir" ]];then
    echo "ERROR: log dir $logdir does not exist"
    exit 1
fi

if [[ ! -d "$tempdir" ]];then
    echo "ERROR: temp dir $tempdir does not exist"
    exit 1
fi

if [[ ! -d "$output" ]];then
    echo "ERROR: root output dir $output does not exist"
    exit 1
fi

# output dir for the current chr
output2="$output"/chr"$c"
mkdir -p "$output2"
if [[ ! -d "$output2" ]];then
    echo "ERROR: could not create output directory $output2 for chr $c"
    exit 1
fi

echo "OUTPUT DIR $output2"
echo "LOG DIR $logdir"
echo "START" $(date)
echo "---------------------------------------------------------------"
echo ""

# file with list of files
flist="$tempdir"/file_list_chr_"$c".txt
echo "" > "$flist"

if [[ "$resume" == "no" ]];then
    # pattern: _c9_ for chrom 9 etc.
    for f in $(find "$input" -maxdepth 1 -mindepth 1 -name "*_c${c}_*.vcf.gz");do
	b=$(basename $f)
	ln -s $f "${output2}/$b"
    done
fi

for f in $(find "${output2}" -mindepth 1 -maxdepth 1 -name "*.vcf.gz");do
    echo "$f" >> "$flist"
done

total=$(cat $flist | wc -l)

if [[ "$total" -eq 0 ]];then
    echo "INFO: no VCFs for chr $c found in $input; exit"
    exit 0
fi

ID=$(sbatch --job-name=collect_stats_chr_"$c" --cpus-per-task=1 --mem-per-cpu=10G --time=10:00:00 -p normal_q --array=1-"$total" -o "$logdir"/collect_stats_chr_"$c"_%A_part_%a.log -e "$logdir"/collect_stats_chr_"$c"_%A_part_%a.err --wrap="singularity exec -B /compute/Genomics /compute/Genomics/containers/worker_3.1 /compute/Genomics/software/scripts/general_purpose_scripts/collect_stats_chunk.sh -i $output2 -p $pheno -f $flist" | cut -d ' ' -f 4)

# remove file list after chr processing is done

sbatch --job-name=cleanup_chr_"$c" --dependency=afterany:"$ID"  --cpus-per-task=1 --mem-per-cpu=1M -p normal_q --array=1 -o "$logdir"/cleanup_chr_"$c"_%A.log -e "$logdir"/cleanup_chr_"$c"_%A.err --wrap="rm -v $flist"
