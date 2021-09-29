#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
script="${scriptdir}/ukbb_select.py"
hesin_script="${scriptdir}/hesin_select.py"
fullsname=$(realpath $0)

function usage () {
    echo ""
    echo "Wrapper script for selecting self reported OA hip/knee/both replacement cases from a UKBB release"
    echo ""
    echo "Usage: OA_select_SR_replacement.sh -r | --release <release of the MAIN dataset>"
    echo "                                   -e | --hesin <release of the HESIN dataset>"
    echo "                                   -o | --output <output prefix>"
    echo "                                   -k | --key <one of: TKR/THR/TJR>"
    echo "                                   -c | --config <optional: config file; default: config.txt in script directory>"
    echo "                                   -h | --help"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

OPTS=$(getopt -o hc:r:o:k:e: -l help,release:,hesin:,config:,output:,key: -n 'OA_select_SR_replacement' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

declare -A keys=( [THR]=1 [TKR]=1 [TJR]=1 )
declare -A keys2=( [THR]=1318 [TKR]=1319 )

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
icd_exclusion_file=""
readValue "$config" OA_CASE_EXCLUSION icd_exclusion_file
exitIfNotFile "$icd_exclusion_file" "ERROR: OA_CASE_EXCLUSION ($icd_exclusion_file) does not exist"

# --------------------------------------------------------------------------------------------------------

logfile="$outprefix".log
outfile="$outprefix".txt
exitIfExists "$outfile" "ERROR: output file $outfile already exists"
: > "$logfile"
out_dir=$(dirname "$outfile")

tmp_ukbb_out=$(mktemp -p "$out_dir" temp_ukbb_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temp UKBB output file"|tee -a "$logfile"
    exit 1
fi
tmp_hesin_out=$(mktemp -p "$out_dir" temp_hesin_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temp HESIN output file"|tee -a "$logfile"
    rm -f "$tmp_ukbb_out"
    exit 1
fi
tmp_icd9=$(mktemp -p "$out_dir" temp_icd9_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temp ICD9 file"|tee -a "$logfile"
    rm -f "$tmp_ukbb_out" "$tmp_hesin_out"
    exit 1
fi
tmp_icd10=$(mktemp -p "$out_dir" temp_icd10_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temp ICD10 file"|tee -a "$logfile"
    rm -f "$tmp_ukbb_out" "$tmp_hesin_out" "$tmp_icd9"
    exit 1
fi

echo "Wrapper script: current dir: ${PWD}" >> "$logfile"
echo "Wrapper script: command line: $scriptname ${args[@]}" >> "$logfile"
echo "" >> "$logfile"

# --------------------------------------------------------------------------------------------------------

PYTHONPATH="${scriptdir}"/python "${script}" -p OA -r "$release" -o "$tmp_ukbb_out" --cc 20002:1465 --cc 20004:1318 --cc 20004:1319 2>>"$logfile"
# select cases based on key

grep ^icd9 "$icd_exclusion_file"| cut -f 2 > "$tmp_icd9"
grep ^icd10 "$icd_exclusion_file"| cut -f 2 > "$tmp_icd10"
PYTHONPATH="${scriptdir}"/python "${hesin_script}" -p OA -r "$hesin_release" -o "$tmp_hesin_out" --icd9 "$tmp_icd9" --icd10 "$tmp_icd10" 2>>"$logfile"

join -1 1 -2 1 -a 1 -t$'\t' -e NULL -o 1.1,2.1 <(tail -n +2 "$tmp_ukbb_out"|grep 1$|cut -f 1|sort) <(sort "$tmp_hesin_out")|grep NULL|cut -f 1 >"$outfile"

rm -f "$tmp_icd9" "$tmp_icd10" "$tmp_ukbb_out" "$tmp_hesin_out"


if [[ "$key" == "TJR" ]];then
    if [[ -z "${outfile}" ]];then
	cat <("$fullsname" -c "$config" -r "$release" -k "TKR" 2>/dev/null) <("$fullsname" -c "$config" -r "$release" -k "THR" 2>/dev/null)|sort|uniq
    else
	cat <("$fullsname" -c "$config" -r "$release" -k "TKR" 2>/dev/null) <("$fullsname" -c "$config" -r "$release" -k "THR" 2>/dev/null)|sort|uniq|gzip - -c >"$outfile" 
    fi
else
    if [[ -z "${outfile}" ]];then 
	"$script" -p "OA" -r "$release" --cc 20002,1465 --cc 20004,"${keys2[$key]}" -c "$config" | tail -n +2 | awk 'BEGIN{FS="\t";}$2=="1" && $3=="1"{print $1;}'
    else
	"$script" -p "OA" -r "$release" --cc 20002,1465 --cc 20004,"${keys2[$key]}" -c "$config" | tail -n +2 | awk 'BEGIN{FS="\t";}$2=="1" && $3=="1"{print $1;}'| gzip - -c >"$outfile" 
    fi
fi

exit 0
