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
        "m" ) mode="${OPTARG}";;
        "t" ) pt="${OPTARG}";;
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
echo "MODE $mode" | ts
echo "PVALUE THRESHOLD $pt" | ts
echo "CURRENT FILENO $n" | ts
echo "CURRENT FILE $fname" | ts
echo "OUTPUT FILE $outname" | ts
echo "--------------------------------------------------------------"
echo ""

pheno_name=$(head -n 1 $pheno| tr ' ' '\t' |cut -f 2)
echo "INFO: phenotype name: $pheno_name"  | ts

if [[ ! -s "$outname" ]];then
    bcftools norm -m- "$fname" | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Ov | qctool2 -g - -filetype vcf -differential "$pheno_name" -osnp "$outname" -s "$pheno"
    retval=$?
else
    echo "INFO: $outname already exists" | ts
fi

if [[ $retval -ne 0 ]];then
    echo "INFO: something went wrong when creating qctool output; exit" | ts
    exit 1
fi

if [[ "$mode" == "stats" ]];then
    if [[ $retval -eq 0 ]];then
	echo "INFO: removing $fname" | ts
	rm -v "$fname"
    else
	echo "INFO: something went wrong; keeping $fname" | ts
    fi
else
    to_remove=${outname/%.qctool.out/.rm}
    # filtering; P-value: lrt_pvalue
    echo "INFO: filtering using P-values" | ts
    echo "INFO: input: $outname" | ts
    echo "INFO: output: $to_remove" | ts
    grep -v "^#" "$outname" | tail -n +2 | awk -v p=${pt} 'BEGIN{FS="\t";OFS="\t";}$13<p{print $2;}' > "$to_remove"
    echo "INFO: done" | ts
    echo "--------------------------------------------------------------"

    final_vcf=${fname/%.vcf.gz/.filtered.vcf.gz}
    # removing variants with P<threshold, merging back
    echo "INFO: bcftools: creating filtered output VCF" | ts
    echo "INFO: input: $fname" | ts
    echo "INFO: input: excluding variants in $to_remove" | ts
    echo "INFO: output: $final_vcf" | ts
    bcftools view --exclude ID=@"$to_remove" "$fname" -Ov | bcftools norm -m+ | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Ov | bcftools +fill-tags -Oz -o "$final_vcf"
    retval=$?
    if [[ $retval -ne 0 ]];then
	echo "INFO: something went wrong when creating filtered VCF; exit" | ts
	exit 1
    fi
    echo "INFO: done" | ts
    echo "--------------------------------------------------------------"

    # TABIX output
    echo "INFO: tabix $final_vcf" | ts
    tabix "$final_vcf"
    if [[ $? -ne 0 ]];then
	echo "INFO: something went wrong when indexing filtered VCF; exit" | ts
	exit 1
    fi
    echo "INFO: done" | ts
    echo "--------------------------------------------------------------"

    echo "INFO: removing intermediate files" | ts
    rm -fv "$fname" "$to_remove" "$outname"
    echo "INFO: done" | ts
fi

exit 0