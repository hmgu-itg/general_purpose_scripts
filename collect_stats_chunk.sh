#!/usr/bin/env bash
#
# using bcftools and qctool2 to collect missingness stats from one VCF file
#
###################

set -eo pipefail

function usage {
    echo ""
    echo "Usage: $0 -i <input dir>"
    echo "          -p <pheno file>"
    echo "          -f <file list>"
    exit 0
}

OPTIND=1
while getopts "i:p:f:" optname; do
    case "$optname" in
        "i" ) input="${OPTARG}";;
        "p" ) pheno="${OPTARG}";;
        "f" ) flist="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

input=${input%/}

# current file number in the file list
n=$SLURM_ARRAY_TASK_ID

fname=$(cat $flist | head -n $n | tail -n 1)
outname=${fname/%.vcf.gz/.qctool.out}

echo "INPUT DIR $input" | ts
echo "PHENOTYPE FILE $pheno" | ts
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
