#!/usr/bin/env bash
#
# using bcftools and qctool2 to collect missingness stats from one VCF file
#
###################

set -eo pipefail

function usage {
    echo ""
    echo "Usage: $0"
    echo "          -p <pheno file>"
    echo "          -f <file list>"
    echo "        { -m : <mode>; optional, \"stats\" or \"full\"; default: \"full\" }"
    echo "        { -t : <pvalue threshold>; only required if mode is \"full\"}"
    exit 0
}

mode="full"
pt=""
OPTIND=1
while getopts "p:f:m:t:" optname; do
    case "$optname" in
        "p" ) pheno="${OPTARG}";;
        "f" ) flist="${OPTARG}";;
        "m" ) input="${OPTARG}";;
        "t" ) input="${OPTARG}";;
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

if [[ ! -f "$flist" ]];then
    echo "ERROR: file list $flist does not exist; exit"
    exit 1
fi

if [[ ! -f "$pheno" ]];then
    echo "ERROR: phenotype file $pheno does not exist; exit"
    exit 1
fi

# current file number in the file list
n=$SLURM_ARRAY_TASK_ID

fname=$(cat $flist | head -n $n | tail -n 1)
outname=${fname/%.vcf.gz/.qctool.out}

echo "INPUT DIR $input" | ts
echo "PHENOTYPE FILE $pheno" | ts
echo "MODE $mode"
echo "PVALUE THRESHOLD $pt"
echo "CURRENT FILENO $n" | ts
echo "CURRENT FILE $fname" | ts
echo "OUTPUT FILE $outname" | ts
echo "--------------------------------------------------------------"
echo ""

pheno_name=$(head -n 1 $pheno| tr ' ' '\t' |cut -f 2)
echo "INFO: phenotype name: $pheno_name"  | ts

bcftools norm -m- "$fname" | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Ov | qctool2 -g - -filetype vcf -differential "$pheno_name" -osnp "$outname" -s "$pheno"

if [[ $? -eq 0 ]];then
    echo "INFO: removing $fname" | ts
    rm -v "$fname"
else
    echo "INFO: something went wrong; keeping $fname" | ts
fi

echo "INFO: done" | ts
echo "--------------------------------------------------------------"
echo ""

exit 0
