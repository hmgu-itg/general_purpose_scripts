#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
hesin_script="${scriptdir}/hesin_select.py"

function usage () {
    echo ""
    echo "Script for selecting hospital diagnosed OA cases from a HESIN release"
    echo ""
    echo "Usage: OA_select_HD.sh -r | --release <release of the HESIN dataset>"
    echo "                       -e | --hesin <release of the HESIN dataset>"
    echo "                       -k | --key <one of: OA, FingerOA, HandOA, HipKneeOA, HipOA, KneeOA, SpineOA, ThumbOA>"
    echo "                       -o | --output <output prefix>"
    echo "                       -c | --config <optional: config file; default: config.txt in script directory>"
    echo "                       -h | --help"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

OPTS=$(getopt -o he:c:r:o:k: -l help,hesin:,release:,config:,output:,key: -n 'OA_select_HD' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

declare -A keys=( [FingerOA]=1 [HandOA]=1 [HipKneeOA]=1 [HipOA]=1 [KneeOA]=1 [OA]=1 [SpineOA]=1 [ThumbOA]=1 )

config=""
release=""
hesin_release=""
outprefix=""
key=""
outfile=""
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
    -r|--release ) release=$2; shift 2 ;;
    -e|--hesin ) hesin_release=$2; shift 2 ;;
    -k|--key ) key=$2; shift 2 ;;
    -c|--config ) config=$2; shift 2 ;;
    -o|--output ) outprefix=$2; shift 2 ;;
    --) shift ; break ;;
    * ) echo "Unexpected option: $1" ; usage ;;
  esac
done

exitIfEmpty "$key" "ERROR: key not specified"
exitIfEmpty "${keys[$key]}" "ERROR: invalid key $key"
exitIfEmpty "$release" "ERROR: MAIN release not specified"
exitIfEmpty "$hesin_release" "ERROR: HESIN release not specified"
exitIfEmpty "$outprefix" "ERROR: output prefix not specified"
if [[ -z "$config" ]];then
    config="${scriptdir}"/config.txt
fi
exitIfNotFile "$config" "ERROR: config $config does not exist"

icd_inclusion_file=""
readValue "$config" HD_OA_ICD icd_inclusion_file
exitIfNotFile "$icd_inclusion_file" "ERROR: HD_OA_ICD ($icd_inclusion_file) does not exist"

icd_exclusion_file=""
readValue "$config" OA_CASE_EXCLUSION icd_exclusion_file
exitIfNotFile "$icd_exclusion_file" "ERROR: OA_CASE_EXCLUSION ($icd_exclusion_file) does not exist"

logfile="$outprefix".log
outfile="$outprefix".txt
exitIfExists "$outfile" "ERROR: output file $outfile already exists"
: > "$logfile"
out_dir=$(dirname "$outfile")

#----------------------------------------------------------------------------------------------------------------

# create temp files
out_dir=$(dirname "$outprefix")
tmp_icd9_excl=$(mktemp -p "$out_dir" temp_icd9_excl_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temp ICD9 file"|tee -a "$logfile"
    exit 1
fi

tmp_icd10_excl=$(mktemp -p "$out_dir" temp_icd10_excl_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temp ICD10 file"|tee -a "$logfile"
    rm -f "$tmp_icd9_excl"
    exit 1
fi

tmp_exclusion_out=$(mktemp -p "$out_dir" temp_icd_excl_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temp ICD exclusion file"|tee -a "$logfile"
    rm -f "$tmp_icd9_excl"
    rm -f "$tmp_icd10_excl"
    exit 1
fi

tmp_icd9_incl=$(mktemp -p "$out_dir" temp_ICD9_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary file"|tee -a "$logfile"
    rm -f "$tmp_icd9_excl"
    rm -f "$tmp_icd10_excl"
    rm -f "$tmp_exclusion_out"
    exit 1
fi

tmp_icd10_incl=$(mktemp -p "$out_dir" temp_ICD10_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary file"|tee -a "$logfile"
    rm -f "$tmp_icd9_incl"
    rm -f "$tmp_icd9_excl"
    rm -f "$tmp_icd10_excl"
    rm -f "$tmp_exclusion_out"
    exit 1
fi

tmp_inclusion_out=$(mktemp -p "$out_dir" temp_inclusion_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary file"|tee -a "$logfile"
    rm -f "$tmp_icd9_incl"
    rm -f "$tmp_icd10_incl"
    rm -f "$tmp_icd9_excl"
    rm -f "$tmp_icd10_excl"
    rm -f "$tmp_exclusion_out"
    exit 1
fi

#----------------------------------------------------------------------------------------------------------------

date "+%d-%b-%Y:%H-%M-%S"|tee -a "$logfile"
echo ""|tee -a "$logfile"
echo "Current dir: ${PWD}"|tee -a "$logfile"
echo "Command line: $scriptname ${args[@]}"|tee -a "$logfile"
echo ""|tee -a "$logfile"
echo "INFO: MAIN release: $release"|tee -a "$logfile"
echo "INFO: HESIN release: $hesin_release"|tee -a "$logfile"
echo "INFO: key: $key"|tee -a "$logfile"
echo "INFO: output prefix: $outprefix"|tee -a "$logfile"
echo "INFO: output file: $outfile"|tee -a "$logfile"
echo ""|tee -a "$logfile"

#----------------------------------------------------------------------------------------------------------------

# exclusion ICD codes
grep ^icd9 "$icd_exclusion_file"| cut -f 2 > "$tmp_icd9_excl"
grep ^icd10 "$icd_exclusion_file"| cut -f 2 > "$tmp_icd10_excl"
PYTHONPATH="${scriptdir}"/python "$hesin_script" -p "OA" -r "$hesin_release" --icd9 "$tmp_icd9_excl" --icd10 "$tmp_icd10_excl" -o "$tmp_exclusion_out" 2>>"$logfile"

# inclusion ICD codes
grep -v "^#" "$icd_inclusion_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd9"{print $3;}' > "$tmp_icd9_incl"
grep -v "^#" "$icd_inclusion_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd10"{print $3;}' > "$tmp_icd10_incl"
PYTHONPATH="${scriptdir}"/python "$hesin_script" -p "OA" -r "$hesin_release" --icd9 "$tmp_icd9_incl" --icd10 "$tmp_icd10_incl" -o "$tmp_inclusion_out" 2>>"$logfile"

# only output samples in inclusion that are not in exclusion
join -1 1 -2 1 -a 1 -t$'\t' -e NULL -o 1.1,2.1 "$tmp_inclusion_out" "$tmp_exclusion_out" | grep NULL | cut -f 1  > "${outfile}"

date "+%d-%b-%Y:%H-%M-%S"|tee -a "$logfile"

# delete temp files
rm -f "$tmp_icd9_incl" "$tmp_icd10_incl" "$tmp_icd9_excl" "$tmp_icd10_excl" "$tmp_exclusion_out" "$tmp_inclusion_out"

exit 0
