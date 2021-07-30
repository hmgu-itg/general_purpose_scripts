#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
selectscript="${scriptdir}/selectCases.pl"

function usage () {
    echo ""
    echo "Script for selecting hospital diagnosed OA cases from a HESIN release"
    echo ""
    echo "Usage: ukbb_select.sh -r | --release <release of the HESIN dataset>"
    echo "                      -k | --key <one of: FingerOA, HandOA, HipKneeOA, HipOA, KneeOA, OA, SpineOA, ThumbOA>"
    echo "                      -o | --output <output prefix; if not specified, output goes to STDOUT>"
    echo "                      -c | --config <optional: config file; default: config.txt in script directory>"
    echo "                      -h | --help"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

OPTS=$(getopt -o hc:r:o:k: -l help,release:,config:,output:,key: -n 'OA_select_HD' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

declare -A available_projects
declare -A keys=( [FingerOA]=1 [HandOA]=1 [HipKneeOA]=1 [HipOA]=1 [KneeOA]=1 [OA]=1 [SpineOA]=1 [ThumbOA]=1 )

config=""
release=""
outprefix=""
key=""
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
    -r|--release ) release=$2; shift 2 ;;
    -k|--key ) key=$2; shift 2 ;;
    -c|--config ) config=$2; shift 2 ;;
    -o|--output ) outprefix=$2; shift 2 ;;
    --) shift ; break ;;
    * ) echo "Unexpected option: $1" ; usage ;;
  esac
done

exitIfEmpty "$key" "ERROR: key not specified"
exitIfEmpty "${keys[$key]}" "ERROR: invalid key $key"
exitIfEmpty "$release" "ERROR: release not specified"
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
readValue "$config" HD_OA_ICD icd_codes_file
exitIfNotFile "$icd_codes_file" "ERROR: HD_OA_ICD not defined in $config"

readValue "$config" DATA_PATH data_path
infile="$data_path"/"${available_projects[$project]}"/releases/hesin_r"$release".tar.gz
exitIfNotFile "$infile" "ERROR: input file $infile does not exist"

tmp_icd9=$(mktemp temp_ICD9_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary file"
    exit 1
fi

tmp_icd10=$(mktemp temp_ICD10_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary file"
    rm "$tmp_icd9"
    exit 1
fi

#----------------------------------------------------------------------------------------------------------------

if [[ ! -z "$outprefix" ]];then
    logfile="$outprefix".log
    outfile="$outprefix".txt.gz
    exitIfExists "$outfile" "ERROR: output file $outfile already exists"
    : > "$logfile"
    sfx="|tee -a $logfile"
else
    sfx="1>&2"
fi

eval date "+%d-%b-%Y:%H-%M-%S" "$sfx"
eval echo "" "$sfx"
eval echo "Current dir: ${PWD}" "$sfx"
eval echo "Command line: $scriptname ${args[@]}" "$sfx"
eval echo "" "$sfx"

eval echo "INFO: release: $release" "$sfx"
eval echo "INFO: key: $key" "$sfx"
eval echo "INFO: input file: $infile" "$sfx"
eval echo "INFO: output prefix: $outprefix" "$sfx"
eval echo "INFO: output file: $outfile" "$sfx"
eval echo "INFO: ICD codes: $icd_codes_file" "$sfx"
eval echo "" "$sfx"

icd9_cn=$(getTGZColNum "${infile}" "hesin_diag.txt" "diag_icd9")
icd10_cn=$(getTGZColNum "${infile}" "hesin_diag.txt" "diag_icd10")
eid_cn=$(getTGZColNum "${infile}" "hesin_diag.txt" "eid")

eval echo "INFO: eid column: ${eid_cn}" "$sfx"
eval echo "INFO: icd9 column: ${icd9_cn}" "$sfx"
eval echo "INFO: icd10 column: ${icd10_cn}" "$sfx"
eval echo "" "$sfx"

indexes=("${icd9_cn}" "${eid_cn}")
readarray -t sorted < <(for a in "${indexes[@]}"; do echo "$a"; done | sort -n)
new_eid1=$(getArrayIndex "${eid_cn}" "${sorted[@]}")
new_icd9=$(getArrayIndex "${icd9_cn}" "${sorted[@]}")

indexes=("${icd10_cn}" "${eid_cn}")
readarray -t sorted < <(for a in "${indexes[@]}"; do echo "$a"; done | sort -n)
new_eid2=$(getArrayIndex "${eid_cn}" "${sorted[@]}")
new_icd10=$(getArrayIndex "${icd10_cn}" "${sorted[@]}")

eval echo "INFO: new eid column: ${new_eid1}" "$sfx"
eval echo "INFO: new icd9 column: ${new_icd9}" "$sfx"
eval echo "" "$sfx"
eval echo "INFO: new eid column: ${new_eid2}" "$sfx"
eval echo "INFO: new icd10 column: ${new_icd10}" "$sfx"
eval echo "" "$sfx"

#----------------------------------------------------------------------------------------------------------------

grep -v "^#" "$icd_codes_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd9"{print $3;}' > "$tmp_icd9"
grep -v "^#" "$icd_codes_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd10"{print $3;}' > "$tmp_icd10"

if [[ -z "${outfile}" ]];then
    cat <(tar -xzf "${infile}" hesin_diag.txt -O|tail -n +2|cut -f "${eid_cn}","${icd9_cn}"|datamash -s -g "${new_eid1}" collapse "$new_icd9"|parallel --pipe --block 10M -N10000 "$selectscript" "$tmp_icd9" 1) <(tar -xzf "${infile}" hesin_diag.txt -O|tail -n +2|cut -f "${eid_cn}","${icd10_cn}"|datamash -s -g "${new_eid2}" collapse "$new_icd10"|parallel --pipe --block 10M -N10000 "$selectscript" "$tmp_icd10" 1) | sort | uniq
else
    cat <(tar -xzf "${infile}" hesin_diag.txt -O|tail -n +2|cut -f "${eid_cn}","${icd9_cn}"|datamash -s -g "${new_eid1}" collapse "$new_icd9"|parallel --pipe --block 10M -N10000 "$selectscript" "$tmp_icd9" 1) <(tar -xzf "${infile}" hesin_diag.txt -O|tail -n +2|cut -f "${eid_cn}","${icd10_cn}"|datamash -s -g "${new_eid2}" collapse "$new_icd10"|parallel --pipe --block 10M -N10000 "$selectscript" "$tmp_icd10" 1) | sort | uniq | gzip - -c > "${outfile}"
fi

eval date "+%d-%b-%Y:%H-%M-%S" "$sfx"

rm "$tmp_icd9" "$tmp_icd10"

exit 0
