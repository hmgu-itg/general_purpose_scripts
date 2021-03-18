#!/usr/bin/env bash
#
# process one VCF chunk
#
###################

function usage {
    echo ""
    echo "Usage: $0 -i <input chunk VCF>"
    echo "          -o <output chunk VCF>"
    echo "          -t <P-value threshold>"
    echo "          -p <pheno file>"
    exit 0
}

t="NA"
OPTIND=1
while getopts "i:o:t:p:" optname; do
    case "$optname" in
        "i" ) input="${OPTARG}";;
        "o" ) output="${OPTARG}";;
        "t" ) t="${OPTARG}";;
        "p" ) pheno="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

if [[ "$t" == "NA" ]];then
    echo "ERROR: P-value threshold not specified"
    usage
    exit 0
fi

if [[ ! -f "$input" ]];then
    echo "ERROR: input VCF file $input does not exist"
    usage
    exit 0
fi

if [[ ! -f "$pheno" ]];then
    echo "ERROR: phenotype file $pheno does not exist"
    usage
    exit 0
fi

# phenotype name from header
pheno_name=$(head -n 1 $pheno| tr ' ' '\t' |cut -f 2)

bcftools norm -m- "$input" | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Ov | qctool2 -g - -filetype vcf -differential "$pheno_name" -osnp difmis_test -s "$pheno"
