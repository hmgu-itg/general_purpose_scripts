#!/usr/bin/env bash

function usage {
    echo ""
    echo "Usage: $0"
    echo "          -i <input dir>"
    echo "          -o <output dir>"
    echo "          -p <pheno file>"
    echo "        { -m : <mode>; optional, \"stats\" or \"full\"; default: \"full\" }"
    echo "        { -t : <pvalue threshold>; only required if mode is \"full\"}"
    echo "        { -r : if analyses should be reusmed; optional; default: false }"
    echo "        { -b : if a pipe should be used instead of temporary files when computing stats. Prevents rerunning.; default: false }"
    echo "        { -c : Collapse multiallelics in output. If set, this flag will create ALT records such as A,C in the output. By default, two records are generated. }"
    echo "        { -s : file containing list of samples to exclude; optional }"
    echo "        { -e <threads> : # threads to use for bcftools/plink; default: 1 }"
    exit 0
}

pipe="no"
resume="no"
mode="full"
pt=""
threads=1
OPTIND=1
collapse="no"

while getopts "i:o:p:m:t:rbcs:e:" optname; do
    case "$optname" in
        "i" ) input="${OPTARG}";;
        "o" ) output="${OPTARG}";;
        "p" ) pheno="${OPTARG}";;
        "m" ) mode="${OPTARG}";;
        "t" ) pt="${OPTARG}";;
        "r" ) resume="yes";;
        "b" ) pipe="yes";;
        "c" ) collapse="yes";;
        "s" ) exsamp="${OPTARG}";;
        "e" ) threads="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    echo "No arguments provided; exit"
    usage
    exit 1
fi

if [[ "$mode" != "full" && "$mode" != "stats" ]];then
    echo "ERROR: mode (-m) should be either \"full\" or \"stats\"; exit"
    exit 1
fi

if [[ "$mode" == "full" && -z "$pt" ]];then
    echo "ERROR: no pvalue threshold (-t) specified; exit"
    exit 1
fi

if [[ ! -f "$pheno" ]];then
    echo "ERROR: phenotype file $pheno does not exist; exit"
    exit 1
fi

input=${input%/}
output=${output%/}

# current chr
c=$SLURM_ARRAY_TASK_ID

echo "INPUT DIR $input"
echo "ROOT OUTPUT DIR $output"
echo "PHENOTYPE FILE $pheno"
echo "MODE $mode"
echo "PVALUE THRESHOLD $pt"
echo "THREADS $threads"
echo "RESUME $resume"
echo "USE PIPE $pipe"
echo "CURRENT CHR $c"

if [[ ! -z "$exsamp" ]];then
      exsamp_opt="-s $exsamp"
      echo "SAMPLE EXCLUSION LIST $exsamp"
fi

if [[ -z "$input" ]];then
    echo "ERROR: input directory not specified; exit"
    exit 1
fi

if [[ -z "$output" ]];then
    echo "ERROR: root output directory not specified; exit"
    exit 1
fi

# should exist already
logdir="$output"/logs
tempdir="$output"/temp

if [[ ! -d "$logdir" ]];then
    echo "ERROR: log dir $logdir does not exist; exit"
    exit 1
fi

if [[ ! -d "$tempdir" ]];then
    echo "ERROR: temp dir $tempdir does not exist; exit"
    exit 1
fi

if [[ ! -d "$output" ]];then
    echo "ERROR: root output dir $output does not exist; exit"
    exit 1
fi

# output dir for the current chr
output2="$output"/chr"$c"

if [[ -d "${output2}" ]];then
    echo "INFO: output directory for chr $c ${output2} already exists"
    if [[ "$resume" == "no" ]];then
	echo "INFO: removing everything in ${output2}"
	rm -vrf "${output2}"/*
    fi
else
    echo "INFO: creating output directory $output for chr $c"
    mkdir -p "${output2}"
fi

if [[ ! -d "${output2}" ]];then
    echo "ERROR: could not create output directory $output2 for chr $c; exit"
    exit 1
fi

echo "OUTPUT DIR $output2"
echo "LOG DIR $logdir"
echo "START" $(date)
echo "---------------------------------------------------------------"
echo ""

# list of files
flist="$tempdir"/file_list_chr_"$c".txt
: > "$flist"

if [[ "$resume" == "no" ]];then
    # pattern: _c9_ for chrom 9 etc.
    # creates symlinks to the input directory VCFs into the output directory
    for f in $(find "$input" -maxdepth 1 -mindepth 1 -name "*_c${c}_*.vcf.gz");do
	b=$(basename $f)
	ln -s $f "${output2}/$b"
    done
fi

# add all VCFs in output dir to the input file list
# If resume=no and the above was run, those are input symlinks
# If resume=yes, filter_chunk will have REMOVED the input file (here a link), so the following will not add to the file
for f in $(find "${output2}" -mindepth 1 -maxdepth 1 -name "*.vcf.gz" -type l);do
    echo "$f" >> "$flist"
done

total=$(cat $flist | wc -l)

# memory per cpu in GB for merging; total: 40G
mem=$(echo "40/$threads" | bc)
if [[ "$mem" -eq 0 ]];then
    mem=1
fi

echo "INFO: MEM $mem"

if [[ "$total" -eq 0 ]];then
    echo "INFO: no input VCFs for chr $c found in $output2"
    rm -fv "$flist"
    if [[ "$mode" == "full" ]];then
	     sbatch --job-name=merge_chunks_chr"$c" --cpus-per-task=${threads} --mem-per-cpu="${mem}"G --time=72:00:00 -p normal_q --array=1 -o "$logdir"/merge_chunks_chr"$c"_%A.log -e "$logdir"/merge_chunks_chr"$c"_%A.err /compute/Genomics/software/scripts/general_purpose_scripts/merge_chunks_chr.sh "$output2" "$threads" "$c"
    fi
else
    topt=""
    if [[ ! -z "$pt" ]];then
	topt="-t $pt"
    fi

    ID=$(sbatch --job-name=filter_chr_"$c" --cpus-per-task=1 --mem-per-cpu=10G --time=10:00:00 -p normal_q --array=1-"$total" -o "$logdir"/filter_chr_"$c"_%A_part_%a.log -e "$logdir"/filter_chr_"$c"_%A_part_%a.err --wrap="singularity exec -B /compute/Genomics /compute/Genomics/containers/worker_3.1 /compute/Genomics/software/scripts/general_purpose_scripts/filter_chunk.sh -p $pheno -f $flist -m $mode $topt" | cut -d ' ' -f 4)

    # remove file list after chr processing is done
    sbatch --job-name=cleanup_chr_"$c" --dependency=afterany:"$ID"  --cpus-per-task=1 --mem-per-cpu=1M -p normal_q --array=1 -o "$logdir"/cleanup_chr_"$c"_%A.log -e "$logdir"/cleanup_chr_"$c"_%A.err --wrap="rm -v $flist"

    if [[ "$mode" == "full" ]];then
	sbatch --job-name=merge_chunks_chr"$c" --dependency=afterok:$ID --cpus-per-task=${threads} --mem-per-cpu="${mem}"G --time=72:00:00 -p normal_q --array=1 -o "$logdir"/merge_chunks_chr"$c"_%A.log -e "$logdir"/merge_chunks_chr"$c"_%A.err /compute/Genomics/software/scripts/general_purpose_scripts/merge_chunks_chr.sh "$output2" "$threads" "$c"
    fi
fi
