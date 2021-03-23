#!/usr/bin/env bash
#
# top level script for processing VCF chunk files
#
########################

function usage {
    echo ""
    echo "Usage: $0"
    echo "          -i <input dir>"
    echo "          -o <output dir>"
    echo "          -p <pheno file>"
    echo "        { -c <chromosome(s)> : optional; default: 1-22 }"
    echo "        { -m : <mode>; optional, \"stats\" or \"full\"; default: \"full\" }"
    echo "        { -t : <pvalue threshold>; only required if mode is \"full\"}"
    echo "        { -r : if analyses should be rerun; optional; default: false }"
    echo "        { -e <threads> : # threads to use for bcftools/plink; default: 1 }"
    exit 0
}

resume="no"
chroms="1-22"
mode="full"
pt=""
threads=1
OPTIND=1
while getopts "i:o:p:c:m:t:re:" optname; do
    case "$optname" in
        "i" ) input="${OPTARG}";;
        "o" ) output="${OPTARG}";;
        "p" ) pheno="${OPTARG}";;
        "c" ) chroms="${OPTARG}";;
        "m" ) mode="${OPTARG}";;
        "t" ) pt="${OPTARG}";;
        "r" ) resume="yes";;
        "e" ) threads="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    echo "No arguments provided; exit"
    usage
    exit 0
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

echo "INPUT DIR $input"
echo "OUTPUT DIR $output"
echo "PHENOTYPE FILE $pheno"
echo "MODE $mode"
echo "PVALUE THRESHOLD $pt"
echo "THREADS $threads"
echo "CHROMS $chroms"
echo "RESUME $resume"
echo "--------------------------------------"
echo ""

logdir="$output"/logs
mkdir -p "$logdir"
if [[ ! -d "$logdir" ]];then
    echo "ERROR: could not create log directory $logdir; exit"
    exit 1
fi

tempdir="$output"/temp
mkdir -p "$tempdir"
if [[ ! -d "$tempdir" ]];then
    echo "ERROR: could not create temporary directory $tempdir; exit"
    exit 1
fi

resopt=""
if [[ "$resume" == "yes" ]];then
    resopt="-r"
fi

topt=""
if [[ ! -z "$pt" ]];then
    topt="-t $pt"
fi

sbatch --job-name=collect_stats --cpus-per-task=1 --mem-per-cpu=1G --time=10:00:00 -p normal_q --array="$chroms" -o "$logdir"/collect_stats_%A_chr_%a.log -e "$logdir"/collect_stats_%A_chr_%a.err /compute/Genomics/software/scripts/general_purpose_scripts/filter_chr.sh -i "$input" -o "$output" -p "$pheno" -m "$mode" -e "$threads" "$topt" "$resopt"
