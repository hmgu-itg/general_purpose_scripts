#!/usr/bin/env bash
#
# top level script for processing VCF chunk files
#
########################

## Get path of current script
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"

if [[ ! -s "$DIR/filter_chr.sh" ]]; then
  echo ERROR: Make sure helper script exists: $DIR/filter_chr.sh
fi

function usage {
    echo ""
    echo "Usage: $0"
    echo "          -i <input dir>"
    echo "          -o <output dir>"
    echo "          -p <pheno file>"
    echo "        { -c <chromosome(s)> : optional; anything that \"sbatch --array=...\" accepts; default: 1-22 }"
    echo "        { -m : <mode>; optional, \"stats\" or \"full\"; default: \"full\" }"
    echo "        { -t : <pvalue threshold>; only required if mode is \"full\"}"
    echo "        { -r : if analyses should be rerun; optional; default: false }"
    echo "        { -b : if a pipe should be used instead of temporary files when computing stats. Prevents rerunning.; default: false }"
    echo "        { -s : file containing list of samples to exclude; optional }"
    echo "        { -e <threads> : optional; # threads to use for bcftools/plink; only required if mode is \"full\"; default: 1 }"
    exit 0
}

resume="no"
pipe="no"
chroms="1-22"
mode="full"
pt=""
threads=1
OPTIND=1
while getopts "i:o:p:c:m:t:rbs:e:" optname; do
    case "$optname" in
        "i" ) input="${OPTARG}";;
        "o" ) output="${OPTARG}";;
        "p" ) pheno="${OPTARG}";;
        "c" ) chroms="${OPTARG}";;
        "m" ) mode="${OPTARG}";;
        "t" ) pt="${OPTARG}";;
        "r" ) resume="yes";;
        "b" ) pipe="yes";;
        "s" ) exsamp="${OPTARG}";;
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
echo "USE PIPE $pipe"
echo "SAMPLE EXCLUSION LIST $exsamp"
echo "--------------------------------------"
echo ""

if [[ -z "$input" ]];then
    echo "ERROR: input directory not specified; exit"
    exit 1
fi

if [[ -z "$output" ]];then
    echo "ERROR: output directory not specified; exit"
    exit 1
fi

if [[ -d "$output" ]];then
    echo "INFO: output directory $output already exists"
else
    echo "INFO: creating output directory $output"
    mkdir -p "$output"
fi

if [[ ! -d "$output" ]];then
    echo "ERROR: could not create output directory $output; exit"
    exit 1
fi

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

pipopt=""
if [[ "$pipe" == "yes" ]];then
    pipopt="-b"
fi


topt=""
if [[ ! -z "$pt" ]];then
    topt="-t $pt"
fi

if [[ ! -z "$exsamp" ]];then
  if [[ ! -s "$exsamp" ]]; then
    echo ERROR: Sample exclusion file does not exist: $exsamp
  fi
      exsamp_opt="-s $exsamp"
fi

sbatch --job-name=filter_chr -c 1 --mem=1G --time=72:00:00 -p normal_q --array="$chroms" -o "$logdir"/filter_chr_%A_chr_%a.log -e "$logdir"/filter_chr_%A_chr_%a.err $DIR/filter_chr.sh -i "$input" -o "$output" -p "$pheno" -m "$mode" -e "$threads" "$topt" "$resopt" "$pipopt" "$exsamp_opt"
