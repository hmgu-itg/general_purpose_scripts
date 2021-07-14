#!/usr/bin/env bash

#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
script="${scriptdir}/selectCases.pl"

function usage () {
    echo ""
    echo "Script for selecting hospital diagnosed OA cases from a HESIN release"
    echo ""
    echo "Usage: ukbb_select.sh -r | --release <release>"
    echo "                      -o | --output <output prefix>"
    echo "                      -c | --config <optional: config file; default: config.txt in script directory>"
    echo "                      -h | --help"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

OPTS=$(getopt -o hc:r:o: -l help,release:,config:,output: -n 'OA_select_HD' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

declare -A available_projects

project="OA"
config=""
release=""
outprefix=""
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
    -r|--release ) release=$2; shift 2 ;;
    -c|--config ) config=$2; shift 2 ;;
    -o|--output ) outprefix=$2; shift 2 ;;
    --) shift ; break ;;
    * ) echo "Unexpected option: $1" ; usage ;;
  esac
done

exitIfEmpty "$release" "ERROR: release not specified"
exitIfEmpty "$outprefix" "ERROR: output prefix not specified"
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

icd9_cn=$(getTGZColNum "${infile}" "hesin_diag.txt" "diag_icd9")
icd10_cn=$(getTGZColNum "${infile}" "hesin_diag.txt" "diag_icd10")
eid_cn=$(getTGZColNum "${infile}" "hesin_diag.txt" "eid")
ins_index_cn=$(getTGZColNum "${infile}" "hesin_diag.txt" "ins_index")

indexes=("${oper4_cn2}" "${eid_cn2}" "${ins_index_cn2}")
readarray -t sorted < <(for a in "${indexes[@]}"; do echo "$a"; done | sort -n)
printf "# %s\n" "${sorted[@]}"
echo "#"

new_eid=$(getArrayIndex "${eid_cn2}" "${sorted[@]}")
new_index=$(getArrayIndex "${ins_index_cn2}" "${sorted[@]}")
new_oper4=$(getArrayIndex "${oper4_cn2}" "${sorted[@]}")

tar -xzf "${infile}" hesin_diag.txt -O|tail -n +2|cut -f "${eid_cn2}","${ins_index_cn2}","${oper4_cn2}"|datamash -s -g "$new_eid","$new_index" collapse "$new_oper4"
