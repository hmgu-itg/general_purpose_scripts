#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
selectscript="${scriptdir}/selectCases.pl"
selectOAHDscript="${scriptdir}/OA_select_HD.sh"
fullsname=$(realpath $0)

function usage () {
    echo ""
    echo "Script for selecting hospital diagnosed OA knee/hip replacement cases from a HESIN release"
    echo ""
    echo "Usage: OA_select_HD_replacement.sh -r | --release <release of the HESIN dataset>"
    echo "                                   -k | --key <one of: THR, TKR, TJR>"
    echo "                                   -o | --output <output prefix>"
    echo "                                   -c | --config <optional: config file; default: config.txt in script directory>"
    echo "                                   -h | --help"
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
declare -A keys2=( [THR]=HipOA [TKR]=KneeOA )

config=""
release=""
outprefix=""
outfile=""
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
exitIfNotFile "$icd_codes_file" "ERROR: HD_OA_ICD does not exist"
readValue "$config" HD_OA_OPCS4 opcs4_codes_file
exitIfNotFile "$opcs4_codes_file" "ERROR: HD_OA_OPCS4 does not exist"

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
eval echo "INFO: OPCS4 codes: $opcs4_codes_file" "$sfx"
eval echo "" "$sfx"

icd9_cn=$(getTGZColNum "${infile}" "hesin_diag.txt" "diag_icd9")
icd10_cn=$(getTGZColNum "${infile}" "hesin_diag.txt" "diag_icd10")
eid_cn=$(getTGZColNum "${infile}" "hesin_diag.txt" "eid")
opcs4_cn=$(getTGZColNum "${infile}" "hesin_oper.txt" "oper4")
index_cn=$(getTGZColNum "${infile}" "hesin_oper.txt" "ins_index")
eid2_cn=$(getTGZColNum "${infile}" "hesin_oper.txt" "eid")

eval echo "INFO: eid \(diag\) column: ${eid_cn}" "$sfx"
eval echo "INFO: icd9 column: ${icd9_cn}" "$sfx"
eval echo "INFO: icd10 column: ${icd10_cn}" "$sfx"
eval echo "INFO: eid \(oper\) column: ${eid2_cn}" "$sfx"
eval echo "INFO: opcs4 column: ${opcs4_cn}" "$sfx"
eval echo "INFO: ins_index \(oper\) column: ${index_cn}" "$sfx"
eval echo "" "$sfx"

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

eval echo "INFO: new eid column: ${new_eid1}" "$sfx"
eval echo "INFO: new icd9 column: ${new_icd9}" "$sfx"
eval echo "" "$sfx"
eval echo "INFO: new eid column: ${new_eid2}" "$sfx"
eval echo "INFO: new icd10 column: ${new_icd10}" "$sfx"
eval echo "" "$sfx"
eval echo "INFO: new eid column \(oper\): ${new_eid_oper}" "$sfx"
eval echo "INFO: new opcs4 column: ${new_opcs4}" "$sfx"
eval echo "INFO: new index column: ${new_index}" "$sfx"
eval echo "" "$sfx"

#----------------------------------------------------------------------------------------------------------------

if [[ "$key" == "TJR" ]];then
    # call itself with TKR and THR and combine (union) the output from both
    if [[ -z "${outfile}" ]];then    
	cat <("$fullsname" -k TKR -r "$release" 2>/dev/null) <("$fullsname" -k THR -r "$release" 2>/dev/null)|sort|uniq
    else
	cat <("$fullsname" -k TKR -r "$release" 2>/dev/null) <("$fullsname" -k THR -r "$release" 2>/dev/null)|sort|uniq|gzip - -c > "${outfile}"
    fi
else
    tmp_opcs4=$(mktemp temp_OPCS4_XXXXXXXX)
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not create temporary file"
	exit 1
    fi
    grep -v "^#" "$opcs4_codes_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k{print $2;}' > "$tmp_opcs4"
    if [[ -z "${outfile}" ]];then
	cat <("$selectOAHDscript" -r "$release" -k "${keys2[$key]}" 2>/dev/null) <(tar -xzf "${infile}" hesin_oper.txt -O|tail -n +2|cut -f "${eid2_cn}","${opcs4_cn}","${index_cn}"|datamash -s -g "${new_eid_oper}","${new_index}" collapse "$new_opcs4"|parallel --pipe --block 10M -N10000 "$selectscript" "$tmp_opcs4" 2|sort|uniq)|sort|uniq -d
    else
	cat <("$selectOAHDscript" -r "$release" -k "${keys2[$key]}" 2>/dev/null) <(tar -xzf "${infile}" hesin_oper.txt -O|tail -n +2|cut -f "${eid2_cn}","${opcs4_cn}","${index_cn}"|datamash -s -g "${new_eid_oper}","${new_index}" collapse "$new_opcs4"|parallel --pipe --block 10M -N10000 "$selectscript" "$tmp_opcs4" 2|sort|uniq)|sort|uniq -d|gzip - -c > "${outfile}"
    fi
    rm "$tmp_opcs4"
fi
    
eval date "+%d-%b-%Y:%H-%M-%S" "$sfx"
exit 0
