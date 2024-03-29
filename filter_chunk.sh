#!/usr/bin/env bash
#
# using bcftools and qctool2 to collect missingness stats from one VCF file
#
###################

set -o pipefail

function usage {
    echo ""
    echo "Usage: $0"
    echo "          -p <pheno file>; Space or tab-separated. Header line 1: ID Phenotype_name; Header line 2: 0 B (for binary)."
    echo "          -f <file list>"
    echo "        { -m : <mode>; optional, \"stats\" or \"full\"; default: \"full\" }"
    echo "        { -t : <pvalue threshold>; only required if mode is \"full\"}"
    echo "        { -o : <output directory>; If not specified, will do scary things with your input files.}"
    echo "        { -n : <chunk number>; If specified, will process the n-th line of the input file list. If unset, will default to SLURM_ARRAY_TASK_ID.}"
    echo "        { -e : <sample exclusion list>; Exclude samples in this file on-the-fly.}"
    echo "        { -b : if a pipe should be used instead of temporary files when computing stats. Prevents rerunning. }"
    echo "        { -c : Collapse multiallelics in output. If set, this flag will create ALT records such as A,C in the output. By default, two records are generated. }"
    echo "        { -x : Just generate exclusion lists and exit. }"
    exit 0
}

mode="full"
pt=""
odir=""
OPTIND=1
pipe="no"
collapse="no"
n=""
elist="no"
justexclude="no"

while getopts "p:f:m:t:o:n:e:bcx" optname; do
    case "$optname" in
        "p" ) pheno="${OPTARG}";;
        "f" ) flist="${OPTARG}";;
        "m" ) mode="${OPTARG}";;
        "t" ) pt="${OPTARG}";;
        "o" ) odir="${OPTARG}";;
        "n" ) n="${OPTARG}";;
        "e" ) elist="${OPTARG}";;
        "b" ) pipe="yes";;
        "c" ) collapse="yes";;
        "x" ) justexclude="yes";;
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

if [[ "$elist" != "no" && ! -f "$flist" ]];then
    echo "ERROR: exclusion list $elist does not exist; exit"
    exit 1
fi


if [[ ! -f "$pheno" ]];then
    echo "ERROR: phenotype file $pheno does not exist; exit"
    exit 1
fi

# current file number in the file list
if [[ -z "$n" ]]; then
  if [[ -z "$SLURM_ARRAY_TASK_ID" ]]; then
    echo "ERROR: no chunk number specified. Use -n or set SLURM_ARRAY_TASK_ID env variable."
    exit 1
  else
    n=$SLURM_ARRAY_TASK_ID
  fi
fi

pheno_name=$(head -n 1 $pheno| tr ' ' '\t' |cut -f 2)
fname=$(cat $flist | head -n $n | tail -n 1)

if [[ -z "$odir" ]]; then
  outname=${fname/%.vcf.gz/."$pheno_name".qctool.out}
else
  obase=$(basename $fname)
  outname=$odir/${obase/%.vcf.gz/.$pheno_name.qctool.out}
fi

if [[ -z "$odir" ]]; then
  opheno=$obase.$pheno.matched
else
  ophenobase=$(basename $pheno)
  opheno=$odir/$obase.$ophenobase.matched
fi


echo "INPUT DIR $input" | ts
echo "PHENOTYPE FILE $pheno" | ts
echo "OUTPUT PHENOTYPE FILE $opheno" | ts

echo "MODE $mode" | ts
echo "PVALUE THRESHOLD $pt" | ts
echo "CURRENT FILENO $n" | ts
echo "CURRENT FILE $fname" | ts
echo "QCTOOL2 OUTPUT FILE $outname" | ts
echo "COLLAPSING VARIANTS $collapse" | ts
echo "PIPE, NO INTERMEDIATE FILES $pipe" | ts
echo "SAMPLE EXCLUSION FILE $efile" | ts
echo "PRODUCE XCL LIST AND EXIT: $justexclude"  | ts
echo "INFO: phenotype name: $pheno_name"  | ts
dmx_vcf=${outname/%.$pheno_name.qctool.out/.dmx.vcf.gz}
echo "INTERMEDIATE VCF $dmx_vcf"
echo "--------------------------------------------------------------"
echo ""

opt_remove_samples=$([ "$elist" == "no" ] && echo "" || echo "$elist")
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

$SCRIPT_DIR/synchronise_VCF_and_sample_file.R $fname $pheno $opheno $opt_remove_samples
retval=$?
if [[ $retval -ne 0 ]];then
    echo "INFO: something went wrong when synchronising phenotype and VCF; exit" | ts
    exit 1
fi

collapse_option() {
  >&2 echo collapse_option called, collapse=$collapse
   if [[ "$collapse" == "yes" ]]; then
       bcftools norm -m +any | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ov
   else
       cat
   fi
 }

 remove_samples_option() {
   >&2 echo "remove_samples_option called, exclusion list set to $elist"
   if [[ "$elist" != "no" ]]; then
     bcftools view -S ^$elist
   else
     cat
   fi
 }

pipe_or_file_input() {
  >&2 echo "pipe_or_file_input called, pipe=$pipe, fname=$fname, dmx_vcf=$dmx_vcf, to_remove=$to_remove"
  if [[ "$pipe" == "yes" ]]; then
    bcftools norm -m- "$fname" -Ov | remove_samples_option | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Oz | bcftools +missing2ref | bcftools view --exclude ID=@"$to_remove"
  else
      bcftools view --exclude ID=@"$to_remove" "$dmx_vcf" -Ov
  fi
}



retval=0

# check if qctool2 output exists
if [[ ! -s "$outname" ]];then
  if [[ "$pipe" == "no" ]]; then
      bcftools norm -m- "$fname" -Ov | remove_samples_option | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ov | bcftools +missing2ref| bgzip > "$dmx_vcf"
      zcat "$dmx_vcf" | qctool2 -g - -filetype vcf -differential "$pheno_name" -osnp "$outname" -s "$opheno"
      retval=$?
    else
      bcftools norm -m- "$fname" -Ov | remove_samples_option | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ov | bcftools +missing2ref | qctool2 -g - -filetype vcf -differential "$pheno_name" -osnp "$outname" -s "$opheno"
      retval=$?
  fi
else
  if [[ "$pipe" == "no" && (! -s "$dmx_vcf") ]]; then
      echo "INFO: $outname already exists, pipe is NO, regenerating $dmx_vcf" | ts
      bcftools norm -m- "$fname" -Ov | remove_samples_option | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ov | bcftools +missing2ref| bgzip > "$dmx_vcf"
  else

    echo "INFO: $outname already exists" | ts
    echo "--------------------------------------------------------------"
    echo ""
  fi
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
    if [[ "$justexclude" == "yes" ]]; then exit 0; fi
    final_vcf=${dmx_vcf/%.dmx.vcf.gz/.filtered."$pheno_name"."$pt".vcf.gz}
    # removing variants with P<threshold, merging back
    echo "INFO: bcftools: creating filtered output VCF" | ts
    echo "INFO: input: $dmx_vcf" | ts
    echo "INFO: input: excluding variants in $to_remove" | ts
    echo "INFO: output: $final_vcf" | ts
    # NORM ANY

    pipe_or_file_input | collapse_option | bcftools +fill-tags -Oz -o "$final_vcf"
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
