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

indir=$(dirname "$input")
if [[ -z "$output" ]];then
    output="$indir"
else
    mkdir -p "$output"
    if [[ ! -d "$output" ]];then
	echo "ERROR: output directory $output could not be created" | ts
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
echo "--------------------------------------------------------------"

x=$(basename $input)

# intermediate files
oname1=${x/%.vcf.gz/.dmx.vcf.gz}
oname1="$output"/"$oname1"
oname2=${x/%.vcf.gz/.qctool.out}
oname2="$output"/"$oname2"
oname3=${x/%.vcf.gz/.Pfilter.out}
oname3="$output"/"$oname3"

# final output
oname4=${x/%.vcf.gz/.filtered.vcf.gz}
oname4="$output"/"$oname4"

# if [[ -s "$oname1" ]];then
#     echo "INFO: bcftools output ($oname1) already exists; skipping" | ts
# else
echo "INFO: demultiplexing using bcftools; input: $input; output: $oname1" | ts
echo "INFO: input: $input" | ts
echo "INFO: output: $oname1" | ts
bcftools norm -m- "$input" | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o "$oname1"
echo "INFO: done" | ts
echo "--------------------------------------------------------------"
echo "INFO: running qctool2; input: $oname1; output: $oname2" | ts
echo "INFO: input: $oname1" | ts
echo "INFO: output: $oname2" | ts
zcat "$oname1" | qctool2 -g - -filetype vcf -differential "$pheno_name" -osnp "$oname2" -s "$pheno"
echo "INFO: done" | ts
echo "--------------------------------------------------------------"
# fi

# if [[ -s "$oname2" ]];then
#     echo "INFO: qctool2 output ($oname2) already exists; skipping" | ts
# else
#     echo "INFO: running qctool2; input: $oname1; output: $oname2" | ts
#     echo "INFO: input: $oname1" | ts
#     echo "INFO: output: $oname2" | ts
#     zcat "$oname1" | qctool2 -g - -filetype vcf -differential "$pheno_name" -osnp "$oname2" -s "$pheno"
#     echo "INFO: done" | ts
#     echo "--------------------------------------------------------------"
# fi

# filtering; P-value: lrt_pvalue
echo "INFO: filtering P-values" | ts
echo "INFO: input: $oname2" | ts
echo "INFO: output: $oname3" | ts
grep -v "^#" "$oname2" | tail -n +2 | awk -v p=$t 'BEGIN{FS="\t";OFS="\t";}$13<p{print $2;}' > "$oname3"
echo "INFO: done" | ts
echo "--------------------------------------------------------------"

# removing variants with P<threshold, merging back
echo "INFO: bcftools: creating filtered output VCF" | ts
echo "INFO: input: $oname1" | ts
echo "INFO: input: excluding variants in $oname3" | ts
echo "INFO: output: $oname4" | ts
bcftools view --exclude ID=@"$oname3" "$oname1" -Ov | bcftools norm -m+ | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Ov | bcftools +fill-tags -Oz -o "$oname4"
echo "INFO: done" | ts
echo "--------------------------------------------------------------"

# TABIX output
echo "INFO: tabix $oname4" | ts
tabix "$oname4"
echo "INFO: done" | ts
echo "--------------------------------------------------------------"

echo "INFO: removing intermediate files" | ts
rm -f "$oname1" "$oname2" "$oname3"
echo "INFO: done" | ts

#rm -rf "$tmp_dir"
