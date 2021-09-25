#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
script="${scriptdir}/ukbb_select.sh"
hesin_script="${scriptdir}/hesin_select.sh"

function usage () {
    echo ""
    echo "Wrapper script for selecting self reported OA cases from a UKBB release"
    echo ""
    echo "Usage: OA_select_SR.sh -r | --release <release of the MAIN dataset>"
    echo "                       -e | --hesin <release of the HESIN dataset>"
    echo "                       -o | --output <output prefix>"
    echo "                       -c | --config <optional: config file; default: config.txt in script directory>"
    echo "                       -h | --help"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

OPTS=$(getopt -o hc:r:o: -l help,release:,config:,output: -n 'OA_select_SR' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

declare -A available_projects
config=""
release=""
hesin_release=""
outprefix=""
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
    -r|--release ) release=$2; shift 2 ;;
    -e|--hesin ) hesin_release=$2; shift 2 ;;
    -c|--config ) config=$2; shift 2 ;;
    -o|--output ) outprefix=$2; shift 2 ;;
    --) shift ; break ;;
    * ) echo "Unexpected option: $1" ; usage ;;
  esac
done

exitIfEmpty "$release" "ERROR: MAIN release not specified"
exitIfEmpty "$hesin_release" "ERROR: HESIN release not specified"
exitIfEmpty "$outprefix" "ERROR: output prefix not specified"
if [[ -z "$config" ]];then
    config="${scriptdir}"/config.txt
fi
exitIfNotFile "$config" "ERROR: config $config does not exist"
readAArray "$config" PROJECTS available_projects
project="OA"
if [[ -z "${available_projects[$project]}" ]];then
    echo "ERROR: project $project is not defined in $config"
    exit 1
fi

icd_exclusion_file=""
readValue "$config" OA_CASE_EXCLUSION icd_exclusion_file
exitIfNotFile "$icd_exclusion_file" "ERROR: OA_CASE_EXCLUSION ($icd_exclusion_file) does not exist"

outfile="${outprefix}".txt
exitIfExists "$outfile" "ERROR: output file $outfile already exists"
logfile="${outprefix}".log
: > "${logfile}"

# create temp files
out_dir=$(dirname "$outfile")
tmp_icd9=$(mktemp -p "$out_dir" temp_icd9_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temp ICD9 file"|tee -a "$logfile"
    exit 1
fi

tmp_icd10=$(mktemp -p "$out_dir" temp_icd10_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temp ICD10 file"|tee -a "$logfile"
    rm -f "$tmp_icd9"
    exit 1
fi

tmp_icd_excl=$(mktemp -p "$out_dir" temp_icd_excl_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temp ICD exclusion file"|tee -a "$logfile"
    rm -f "$tmp_icd9"
    rm -f "$tmp_icd10"
    exit 1
fi

tmp_sr_oa=$(mktemp -p "$out_dir" temp_sr_oa_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temp SR OA inclusion file"|tee -a "$logfile"
    rm -f "$tmp_icd9"
    rm -f "$tmp_icd10"
    rm -f "$tmp_icd_excl"
    exit 1
fi

# save exclusion ICD codes
grep ^icd9 "$icd_exclusion_file"| cut -f 2 > "$tmp_icd9"
grep ^icd10 "$icd_exclusion_file"| cut -f 2 > "$tmp_icd10"

# save IDs having exclusion ICD codes
"$hesin_script" -p "OA" -r "$hesin_release" --icd9 "$tmp_icd9" --icd10 "$tmp_icd10" 2>>"$logfile" | sort >"$tmp_icd_excl" 

# save self-reported OA cases
"$script" -p "OA" -r "$release" --cc 20002,1465 -c "$config" 2>>"$logfile" | tail -n +2 | grep "1$" | cut -f 1 | sort >"$tmp_sr_oa"

# output self-reported OA cases without exclusion ICD codes
join -1 1 -2 1 -a 1 -t$'\t' -e NULL -o 1.1,2.1 "$tmp_sr_oa" "$tmp_icd_excl" | grep NULL | cut -f 1 >"$outfile"

# delete temp files
rm -f "$tmp_icd9" "$tmp_icd10" "$tmp_icd_excl" "$tmp_sr_oa"

exit 0
