#!/usr/bin/env bash
#
# using bcftools and qctool2 to collect missingness stats from one VCF file
#
###################

set -o pipefail

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

pheno_name=$(head -n 1 $pheno| tr ' ' '\t' |cut -f 2)
fname=$(cat $flist | head -n $n | tail -n 1)
outname=${fname/%.vcf.gz/."$pheno_name".qctool.out}

echo "INPUT DIR $input" | ts
echo "PHENOTYPE FILE $pheno" | ts
echo "MODE $mode" | ts
echo "PVALUE THRESHOLD $pt" | ts
echo "CURRENT FILENO $n" | ts
echo "CURRENT FILE $fname" | ts
echo "QCTOOL2 OUTPUT FILE $outname" | ts
echo "--------------------------------------------------------------"
echo ""

echo "INFO: phenotype name: $pheno_name"  | ts
dmx_vcf=${fname/%.vcf.gz/.dmx.vcf.gz}

retval=0

# check if qctool2 output exists
if [[ ! -s "$outname" ]];then
    bcftools norm -m- "$fname" -Ov | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o "$dmx_vcf"
    zcat "$dmx_vcf" | qctool2 -g - -filetype vcf -differential "$pheno_name" -osnp "$outname" -s "$pheno"
    retval=$?
else
    echo "INFO: $outname already exists" | ts
    echo "--------------------------------------------------------------"
    echo ""
fi

if [[ $retval -ne 0 ]];then
    echo "INFO: something went wrong when creating qctool output; exit" | ts
    if [[ -f "$outname" ]];then rm -v "$outname";fi
    if [[ -f "$dmx_vcf" ]];then rm -v "$dmx_vcf";fi
    exit 1
fi

if [[ "$mode" == "stats" ]];then
    echo "INFO: removing $dmx_vcf" | ts
    rm -v "$dmx_vcf"
else
    to_remove=${outname/%.qctool.out/."$pt".rm}
    # filtering; P-value: lrt_pvalue
    echo "INFO: filtering using P-values" | ts
    echo "INFO: input: $outname" | ts
    echo "INFO: output: $to_remove" | ts
    # in case input VCF has incorrect GT records there will be error messages in the qctool output
    echo "GREP"
    grep -i error "$outname" | cut -f 2 > "$to_remove"
    echo "AWK"
    grep -v "^#" "$outname" | tail -n +2 | grep -v -i error | awk -v p=${pt} 'BEGIN{FS="\t";OFS="\t";}$13<p{print $2;}' >> "$to_remove"
    echo "SPONGE"
    sort "$to_remove" | uniq | sponge "$to_remove"
    echo "SPONGE DONE"
    echo "INFO: done" | ts
    echo "--------------------------------------------------------------"

    final_vcf=${fname/%.vcf.gz/.filtered."$pheno_name"."$pt".vcf.gz}
    # removing variants with P<threshold, merging back
    echo "INFO: bcftools: creating filtered output VCF" | ts
    echo "INFO: input: $dmx_vcf" | ts
    echo "INFO: input: excluding variants in $to_remove" | ts
    echo "INFO: output: $final_vcf" | ts
    # NORM ANY
    bcftools view --exclude ID=@"$to_remove" "$dmx_vcf" -Ov | bcftools norm -m +any | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Ov | bcftools +fill-tags -Oz -o "$final_vcf"
    retval=$?
    if [[ $retval -ne 0 ]];then
	echo "INFO: something went wrong when creating filtered VCF; exit" | ts
	if [[ -f "$final_vcf" ]];then rm -v "$final_vcf";fi
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
    rm -fv "$fname" "$to_remove" "$outname" "$dmx_vcf"
fi

echo "INFO: done" | ts

exit 0
