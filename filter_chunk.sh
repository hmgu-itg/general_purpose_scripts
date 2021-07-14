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
    echo "        { -o : <output directory>; If not specified, will do scary things with your input files.}"
    echo "        { -b : if a pipe should be used instead of temporary files when computing stats. Prevents rerunning. }"
    echo "        { -c : Collapse multiallelics in output. If set, this flag will create ALT records such as A,C in the output. By default, two records are generated. }"
    exit 0
}

mode="full"
pt=""
odir=""
OPTIND=1
pipe="no"
collapse="no"

while getopts "p:f:m:t:o:bc" optname; do
    case "$optname" in
        "p" ) pheno="${OPTARG}";;
        "f" ) flist="${OPTARG}";;
        "m" ) mode="${OPTARG}";;
        "t" ) pt="${OPTARG}";;
        "o" ) odir="${OPTARG}";;
        "b" ) pipe="yes";;
        "c" ) collapse="yes";;
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

if [[ ! -z "$exsamp" ]];then
      exsamp_opt="--excl-samples $exsamp"
fi

# current file number in the file list
n=$SLURM_ARRAY_TASK_ID

pheno_name=$(head -n 1 $pheno| tr ' ' '\t' |cut -f 2)
fname=$(cat $flist | head -n $n | tail -n 1)

if [[ -z "$odir" ]]; then
  outname=${fname/%.vcf.gz/."$pheno_name".qctool.out}
else
  obase=$(basename $fname)
  outname=$odir/${obase/%.vcf.gz/.$pheno_name.qctool.out}
fi

echo "INPUT DIR $input" | ts
echo "PHENOTYPE FILE $pheno" | ts
echo "MODE $mode" | ts
echo "PVALUE THRESHOLD $pt" | ts
echo "CURRENT FILENO $n" | ts
echo "CURRENT FILE $fname" | ts
echo "QCTOOL2 OUTPUT FILE $outname" | ts
echo "COLLAPSING VARIANTS $collapse" | ts
echo "PIPE, NO INTERMEDIATE FILES $pipe" | ts
echo "--------------------------------------------------------------"
echo ""

echo "INFO: phenotype name: $pheno_name"  | ts
dmx_vcf=${fname/%.vcf.gz/.dmx.vcf.gz}

collapse_option() {
  echo collapse_option called, collapse=$collapse
   if [[ "$collapse" == "yes" ]]; then
       bcftools norm -m +any | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Ov
   else
       cat
   fi
 }

pipe_or_file_input() {
  echo "pipe_or_file_input called, pipe=$pipe, fname=$fname, dmx_vcf=$dmx_vcf, to_remove=$to_remove"
  if [[ "$pipe" == "yes" ]]; then
    bcftools norm -m- "$fname" -Ov | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o "$dmx_vcf"
  else
      bcftools view --exclude ID=@"$to_remove" "$dmx_vcf" -Ov
  fi
}


retval=0

# check if qctool2 output exists
if [[ ! -s "$outname" ]];then
  if [[ "$pipe" == "no" ]]; then
      bcftools norm -m- "$fname" -Ov | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o "$dmx_vcf"
      zcat "$dmx_vcf" | qctool2 -g - -filetype vcf -differential "$pheno_name" -osnp "$outname" -s "$pheno" $exsamp_opt
      retval=$?
    else
      bcftools norm -m- "$fname" -Ov | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' | qctool2 -g - -filetype vcf -differential "$pheno_name" -osnp "$outname" -s "$pheno" $exsamp_opt
      retval=$?
  fi
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
    echo "INFO: done" | ts
    # echo "INFO: removing $dmx_vcf" | ts
    # rm -v "$dmx_vcf"
else
    to_remove=${outname/%.qctool.out/."$pt".rm}
    # filtering; P-value: lrt_pvalue
    echo "INFO: filtering using P-values" | ts
    echo "INFO: input: $outname" | ts
    echo "INFO: output: $to_remove" | ts
    # in case input VCF has incorrect GT records there will be error messages in the qctool output
    grep -i error "$outname" | cut -f 2 > "$to_remove"
    grep -v "^#" "$outname" | tail -n +2 | grep -v -i error | awk -v p=${pt} 'BEGIN{FS="\t";OFS="\t";}$13<p{print $2;}' >> "$to_remove"
    sort "$to_remove" | uniq | sponge "$to_remove"
    echo "INFO: done" | ts
    echo "--------------------------------------------------------------"

    final_vcf=${fname/%.vcf.gz/.filtered."$pheno_name"."$pt".vcf.gz}
    # removing variants with P<threshold, merging back
    echo "INFO: bcftools: creating filtered output VCF" | ts
    echo "INFO: input: $dmx_vcf" | ts
    echo "INFO: input: excluding variants in $to_remove" | ts
    echo "INFO: output: $final_vcf" | ts
    # NORM ANY

    bcftools view --exclude ID=@"$to_remove" "$dmx_vcf" -Ov | collapse_option() | bcftools +fill-tags -Oz -o "$final_vcf"
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

    ## REMOVED THIS.
    ## This is the last in a line of 3 scripts that actually relies on deleting its input files to work well
    ## (it assumes inputs are all links)

    ## TODO
    ## The whole "resume" section needs to be reworked. No test or linking should be done on the input files
    ## Currently "mode=stats" just generates statistics, and if used with "pipe=no" can be used to "resume" if "mode=full"
    ## Just test the intermediate files for resume: the intermediate VCFs and exclusion lists need to be there. If they're there proceed to filtering and chr merging.

    # echo "INFO: removing intermediate files" | ts
    # rm -fv "$fname" "$to_remove" "$outname" "$dmx_vcf"
    echo "INFO: done" | ts
fi

exit 0
