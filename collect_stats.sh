#!/usr/bin/env bash
#
# top level script for collecting missingness stats
#
########################

function usage {
    echo ""
    echo "Usage: $0"
    echo "          -i <input dir>"
    echo "          -o <output dir>"
    echo "          -p <pheno file>"
    echo "        { -c <chromosome(s)> : optional; default: 1-22 }"
    echo "        { -r : if analyses should be reusmed; optional; default: false }"
    exit 0
}

resume="no"
chroms="1-22"
OPTIND=1
while getopts "i:o:p:c:r" optname; do
    case "$optname" in
        "i" ) input="${OPTARG}";;
        "o" ) output="${OPTARG}";;
        "p" ) pheno="${OPTARG}";;
        "c" ) chroms="${OPTARG}";;
        "r" ) resume="yes";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    echo "No arguments provided; exit"
    usage
    exit 0
fi

input=${input%/}
output=${output%/}

echo "INPUT DIR $input"
echo "OUTPUT DIR $output"
echo "PHENOTYPE FILE $pheno"
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

sbatch --job-name=collect_stats --cpus-per-task=1 --mem-per-cpu=1G --time=10:00:00 -p normal_q --array="$chroms" -o "$logdir"/collect_stats_%A_chr_%a.log -e "$logdir"/collect_stats_%A_chr_%a.err /compute/Genomics/software/scripts/general_purpose_scripts/collect_stats_chr.sh -i "$indir" -o "$outdir" -p "$pheno" "$resopt"
