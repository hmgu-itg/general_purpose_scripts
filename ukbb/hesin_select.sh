#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

function usage () {
    echo ""
    echo "Script for selecting samples from HESIN tables"
    echo "Select samples that have at least one ICD code from the supplied list(s)"
    echo ""
    echo "Usage: hesin_select.sh -p | --project <project name>"
    echo "                       -r | --release <release>"
    echo "                       --icd9 <list of ICD9 codes>"
    echo "                       --icd10 <list of ICD10 codes>"
    echo "                      -o | --output <optional: output prefix; if not specified, output goes to STDOUT>"
    echo "                      -c | --config <optional: config file; default: config.txt in script directory>"
    echo "                      -h | --help"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

OPTS=$(getopt -o hc:p:r:o: -l help,project:,release:,config:,icd9:,icd10:,output: -n 'hesin_select' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

#----------------------------------------------------------------------------------------------------------------

declare -A available_projects

config=""
release=""
icd9_file=""
icd10_file=""
outprefix=""
outfile=""
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
    -r|--release ) release=$2; shift 2 ;;
    -p|--project ) project=$2; shift 2 ;;
    -c|--config ) config=$2; shift 2 ;;
    --icd9 ) icd9_file=$2; shift 2 ;;
    --icd10 ) icd10_file=$2; shift 2 ;;
    -o|--output ) outprefix=$2; shift 2 ;;
    --) shift ; break ;;
    * ) echo "Unexpected option: $1" ; usage ;;
  esac
done

# at least one ICD input is needed

exitIfEmpty "$project" "ERROR: project not specified"
exitIfEmpty "$release" "ERROR: release not specified"
if [[ -z "$config" ]];then
    config="${scriptdir}"/config.txt
fi
exitIfNotFile "$config" "ERROR: config $config does not exist"
readAArray "$config" PROJECTS available_projects
if [[ -z "${available_projects[$project]}" ]];then
    echo "ERROR: project $project is not defined in $config"
    exit 1
fi
readValue "$config" DATA_PATH data_path
infile="$data_path"/"${available_projects[$project]}"/releases/hesin_r"$release".tar.gz
exitIfNotFile "$infile" "ERROR: input file $infile does not exist"

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

#----------------------------------------------------------------------------------------------------------------

if [[ -z "${outfile}" ]];then
    cat <(tar -xzf "${infile}" hesin_diag.txt -O|tail -n +2|cut -f "${eid_cn}","${icd9_cn}"|sort -k2,2|join -t$'\t' -1 2 -2 1 - <(sort ${icd9_file}|uniq)|cut -f 2|sort|uniq) <(tar -xzf "${infile}" hesin_diag.txt -O|tail -n +2|cut -f "${eid_cn}","${icd10_cn}"|sort -k2,2|join -t$'\t' -1 2 -2 1 - <(sort ${icd10_file}|uniq)|cut -f 2|sort|uniq)|sort|uniq
else
    cat <(tar -xzf "${infile}" hesin_diag.txt -O|tail -n +2|cut -f "${eid_cn}","${icd9_cn}"|sort -k2,2|join -t$'\t' -1 2 -2 1 - <(sort ${icd9_file}|uniq)|cut -f 2|sort|uniq) <(tar -xzf "${infile}" hesin_diag.txt -O|tail -n +2|cut -f "${eid_cn}","${icd10_cn}"|sort -k2,2|join -t$'\t' -1 2 -2 1 - <(sort ${icd10_file}|uniq)|cut -f 2|sort|uniq)|sort|uniq|gzip - -c > "${outfile}"
fi

eval date "+%d-%b-%Y:%H-%M-%S" "$sfx"
