#!/usr/bin/env bash
#
# process one VCF chunk
#
###################

function usage {
    echo ""
    echo "Usage: $0 -i <input chunk VCF>"
    echo "          -t <P-value threshold>"
    echo "          -p <pheno file>"
    echo "        { -o <output dir> : optional; default: input file directory }"
    exit 0
}

t="NA"
OPTIND=1
output=""
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
    echo "ERROR: P-value threshold not specified" | ts
    usage
    exit 0
fi

if [[ ! -f "$pheno" ]];then
    echo "ERROR: phenotype file $pheno does not exist" | ts
    usage
    exit 0
fi

if [[ ! -f "$input" ]];then
    echo "ERROR: input VCF file $input does not exist" | ts
    usage
    exit 0
fi

indir=$(dirname $(realpath "$input"))
if [[ -z "$output" ]];then
    output="$indir"
else
    if [[ ! -d "$output" ]];then
	echo "ERROR: output directory $output does not exist" | ts
	usage
	exit 0
    fi
fi

# dname=$(dirname $(realpath "$input"))
# tmp_dir=$(mktemp -d -p "$dname" -t process_chunk-XXXXXXX)
# if [[ ! -d "$tmp_dir" ]];then
#     echo "ERROR: could not create temporary directory $tmp_dir" | ts
#     exit 1
# else
#     echo "INFO: created temporary directory $tmp_dir" | ts
# fi

# phenotype name from header
pheno_name=$(head -n 1 $pheno| tr ' ' '\t' |cut -f 2)
echo "INFO: phenotype name: $pheno_name"  | ts

x=$(basename $input)
oname1=${x/%.vcf.gz/.dmx.vcf.gz}
oname1="$output"/"$oname1"
oname2=${x/%.vcf.gz/.qctool.out}
oname2="$output"/"$oname2"
oname3=${x/%.vcf.gz/.Pfilter.out}
oname3="$output"/"$oname3"

if [[ -s "$oname1" ]];then
    echo "INFO: bcftools output ($oname1) already exists; skipping" | ts
else
    echo "INFO: running bcftools; input: $input; output: $oname1" | ts
    bcftools norm -m- "$input" | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o "$oname1"
fi

if [[ -s "$oname2" ]];then
    echo "INFO: qctool2 output ($oname2) already exists; skipping" | ts
else
    echo "INFO: running qctool2; input: $oname1; output: $oname2" | ts
    zcat "$oname1" | qctool2 -g - -filetype vcf -differential "$pheno_name" -osnp "$oname2" -s "$pheno"
fi

# P-value: lrt_pvalue
echo "INFO: filtering P-values" | ts
cat "$oname2" | grep -v "^#" | tail -n +2 | awk -v p=$t 'BEGIN{FS="\t";OFS="\t";}$13<p{print $2;}' > "$oname3"


#rm -rf "$tmp_dir"
