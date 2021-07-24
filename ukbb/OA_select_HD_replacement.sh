#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
selectscript="${scriptdir}/selectCases.pl"

function usage () {
    echo ""
    echo "Script for selecting hospital diagnosed OA knee/hip replacement cases from a HESIN release"
    echo ""
    echo "Usage: ukbb_select.sh -r | --release <release of the HESIN dataset>"
    echo "                      -k | --key <one of: THR, TKR, TJR>"
    echo "                      -o | --output <output prefix>"
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
declare -A keys=( [THR]=1 [TKR]=1 [TJR]=1 )

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
exitIfEmpty "${keys[$key]}" "ERROR: wrong key $key"
outfile="${outprefix}".txt.gz
logfile="${outprefix}".log
exitIfExists "$outfile" "ERROR: output file $outfile already exists"
exitIfEmpty "$release" "ERROR: release not specified"
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
readValue "$config" HD_OA_ICD icd_codes_file
exitIfNotFile "$icd_codes_file" "ERROR: HD_OA_ICD does not exist"
readValue "$config" HD_OA_OPCS4 opcs4_codes_file
exitIfNotFile "$opcs4_codes_file" "ERROR: HD_OA_OPCS4 does not exist"

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

tmp_opcs4=$(mktemp temp_OPCS4_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary file"
    rm "$tmp_icd9"
    rm "$tmp_icd10"
    exit 1
fi

#----------------------------------------------------------------------------------------------------------------

: > "$logfile"

date "+%d-%b-%Y:%H-%M-%S" | tee -a "$logfile"
echo "" | tee -a "$logfile"
echo "Current dir: ${PWD}" | tee -a "$logfile"
echo "Command line: $scriptname ${args[@]}" | tee -a "$logfile"
echo "" | tee -a "$logfile"

echo "INFO: release: $release" | tee -a "$logfile"
echo "INFO: key: $key" | tee -a "$logfile"
echo "INFO: input file: $infile" | tee -a "$logfile"
echo "INFO: output prefix: $outprefix" | tee -a "$logfile"
echo "INFO: output file: $outfile" | tee -a "$logfile"
echo "INFO: ICD codes: $icd_codes_file" | tee -a "$logfile"
echo "INFO: OPCS4 codes: $opcs4_codes_file" | tee -a "$logfile"
echo "" | tee -a "$logfile"

icd9_cn=$(getTGZColNum "${infile}" "hesin_diag.txt" "diag_icd9")
icd10_cn=$(getTGZColNum "${infile}" "hesin_diag.txt" "diag_icd10")
eid_cn=$(getTGZColNum "${infile}" "hesin_diag.txt" "eid")
opcs4_cn=$(getTGZColNum "${infile}" "hesin_oper.txt" "oper4")
index_cn=$(getTGZColNum "${infile}" "hesin_oper.txt" "ins_index")
eid2_cn=$(getTGZColNum "${infile}" "hesin_oper.txt" "eid")

echo "INFO: eid (diag) column: ${eid_cn}" | tee -a "$logfile"
echo "INFO: icd9 column: ${icd9_cn}" | tee -a "$logfile"
echo "INFO: icd10 column: ${icd10_cn}" | tee -a "$logfile"
echo "INFO: eid (oper) column: ${eid2_cn}" | tee -a "$logfile"
echo "INFO: opcs4 column: ${opcs4_cn}" | tee -a "$logfile"
echo "INFO: ins_index (oper) column: ${index_cn}" | tee -a "$logfile"
echo "" | tee -a "$logfile"

indexes=("${icd9_cn}" "${eid_cn}")
readarray -t sorted < <(for a in "${indexes[@]}"; do echo "$a"; done | sort -n)
new_eid1=$(getArrayIndex "${eid_cn}" "${sorted[@]}")
new_icd9=$(getArrayIndex "${icd9_cn}" "${sorted[@]}")

indexes=("${icd10_cn}" "${eid_cn}")
readarray -t sorted < <(for a in "${indexes[@]}"; do echo "$a"; done | sort -n)
new_eid2=$(getArrayIndex "${eid_cn}" "${sorted[@]}")
new_icd10=$(getArrayIndex "${icd10_cn}" "${sorted[@]}")

indexes=("${opcs4_cn}" "${index_cn}" "${eid2_cn}")
readarray -t sorted < <(for a in "${indexes[@]}"; do echo "$a"; done | sort -n)
new_eid_oper=$(getArrayIndex "${eid2_cn}" "${sorted[@]}")
new_opcs4=$(getArrayIndex "${opcs4_cn}" "${sorted[@]}")
new_index=$(getArrayIndex "${index_cn}" "${sorted[@]}")

echo "INFO: new eid column: ${new_eid1}" | tee -a "$logfile"
echo "INFO: new icd9 column: ${new_icd9}" | tee -a "$logfile"
echo "" | tee -a "$logfile"
echo "INFO: new eid column: ${new_eid2}" | tee -a "$logfile"
echo "INFO: new icd10 column: ${new_icd10}" | tee -a "$logfile"
echo "" | tee -a "$logfile"
echo "INFO: new eid column (oper): ${new_eid_oper}" | tee -a "$logfile"
echo "INFO: new opcs4 column: ${new_opcs4}" | tee -a "$logfile"
echo "INFO: new index column: ${new_index}" | tee -a "$logfile"
echo "" | tee -a "$logfile"

#----------------------------------------------------------------------------------------------------------------

grep -v "^#" "$icd_codes_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd9"{print $3;}' > "$tmp_icd9"
grep -v "^#" "$icd_codes_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd10"{print $3;}' > "$tmp_icd10"
grep -v "^#" "$opcs4_codes_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k{print $2;}' > "$tmp_opcs4"

cat <(cat <(tar -xzf "${infile}" hesin_diag.txt -O|tail -n +2|cut -f "${eid_cn}","${icd9_cn}"|datamash -s -g "${new_eid1}" collapse "$new_icd9"|parallel --pipe --block 10M -N10000 "$selectscript" "$tmp_icd9" 1) <(tar -xzf "${infile}" hesin_diag.txt -O|tail -n +2|cut -f "${eid_cn}","${icd10_cn}"|datamash -s -g "${new_eid2}" collapse "$new_icd10"|parallel --pipe --block 10M -N10000 "$selectscript" "$tmp_icd10" 1)|sort|uniq) <(tar -xzf "${infile}" hesin_oper.txt -O|tail -n +2|cut -f "${eid2_cn}","${opcs4_cn}","${index_cn}"|datamash -s -g "${new_eid_oper}","${new_index}" collapse "$new_opcs4"|parallel --pipe --block 10M -N10000 "$selectscript" "$tmp_opcs4" 2)|sort|uniq -d|gzip - -c > "${outfile}"

date "+%d-%b-%Y:%H-%M-%S" | tee -a "$logfile"
rm "$tmp_icd9" "$tmp_icd10" "$tmp_opcs4"
exit 0
