#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
upperdir=$(dirname $scriptdir)
source "${upperdir}/functions.sh"
script="${upperdir}/ukbb_select.py"
hesin_script="${upperdir}/hesin_select.py"

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
    config="${upperdir}"/config.txt
fi
exitIfNotFile "$config" "ERROR: config $config does not exist"
icd_exclusion_file=""
readValue "$config" OA_CASE_EXCLUSION icd_exclusion_file
exitIfNotFile "$icd_exclusion_file" "ERROR: OA_CASE_EXCLUSION ($icd_exclusion_file) does not exist"

logfile="$outprefix".log
outfile="$outprefix".txt
exitIfExists "$outfile" "ERROR: output file $outfile already exists"
: > "$logfile"
out_dir=$(dirname "$outfile")

# --------------------------------------------------------------------------------------------------------

# create temp files
declare -A tempfiles
tempfiles["tmp_ukbb_out"]="tmp_XXXXXXX"
tempfiles["tmp_hesin_out"]="tmp_XXXXXXX"
tempfiles["tmp_icd9"]="tmp_XXXXXXX"
tempfiles["tmp_icd10"]="tmp_XXXXXXX"

if createTempFiles "$out_dir" tempfiles;then
    echo "INFO: created temporary files"|tee -a "$logfile"
else
    echo "ERROR: failed to create temporary files"|tee -a "$logfile"
    exit 1
fi

# --------------------------------------------------------------------------------------------------------

echo "Wrapper script: current dir: ${PWD}" | tee -a "$logfile"
echo "Wrapper script: command line: $scriptname ${args[@]}" | tee -a "$logfile"
echo "" | tee -a "$logfile"
echo "Wrapper script: MAIN release: $release"  | tee -a "$logfile"
echo "Wrapper script: HESIN release: $hesin_release"  | tee -a "$logfile"
echo "Wrapper script: key: $key"  | tee -a "$logfile"
echo "Wrapper script: output prefix: $outprefix"  | tee -a "$logfile"
echo "Wrapper script: output file: $outfile"  | tee -a "$logfile"
echo ""  | tee -a "$logfile"

# --------------------------------------------------------------------------------------------------------
echo "RUNNING ${script}" | tee -a "$logfile"
echo "" | tee -a "$logfile"

# CASE INCLUSION
# first we create an output file with eid and three case-control columns for SR OA, SR hip replacement and SR knee replacement
PYTHONPATH="${upperdir}"/python "${script}" -p OA -r "$release" -o "${tempfiles[tmp_ukbb_out]}" --cc 20002:1465 --cc 20004:1318 --cc 20004:1319 > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)

# select cases based on the key: THR, TKR or TJR
declare -A colnumbers
getColNumbers "${tempfiles[tmp_ukbb_out]}" "cat" colnumbers
if [[ "$key" == "THR" ]];then
    tail -n +2 "${tempfiles[tmp_ukbb_out]}" | awk -v a="${colnumbers[f.eid]}" -v b="${colnumbers[cc-20002-1465]}" -v c="${colnumbers[cc-20004-1318]}" 'BEGIN{FS="\t";}$b=="1" && $c=="1"{print $a;}' | sponge "${tempfiles[tmp_ukbb_out]}"
fi
if [[ "$key" == "TKR" ]];then
    tail -n +2 "${tempfiles[tmp_ukbb_out]}" | awk -v a="${colnumbers[f.eid]}" -v b="${colnumbers[cc-20002-1465]}" -v c="${colnumbers[cc-20004-1319]}" 'BEGIN{FS="\t";}$b=="1" && $c=="1"{print $a;}' | sponge "${tempfiles[tmp_ukbb_out]}"
fi
if [[ "$key" == "TJR" ]];then
    tail -n +2 "${tempfiles[tmp_ukbb_out]}" | awk -v a="${colnumbers[f.eid]}" -v b="${colnumbers[cc-20002-1465]}" -v c="${colnumbers[cc-20004-1318]}" -v d="${colnumbers[cc-20004-1319]}" 'BEGIN{FS="\t";}$b=="1" && ($c=="1" || $d=="1"){print $a;}'|sponge "${tempfiles[tmp_ukbb_out]}"
fi

n=$(cat "${tempfiles[tmp_ukbb_out]}"| wc -l)
echo "" | tee -a "$logfile"
echo "INFO: $n samples meet inclusion criteria" | tee -a "$logfile"

# CASE EXCLUSION
grep ^icd9 "$icd_exclusion_file"| cut -f 2 > "${tempfiles[tmp_icd9]}"
grep ^icd10 "$icd_exclusion_file"| cut -f 2 > "${tempfiles[tmp_icd10]}"

echo "" | tee -a "$logfile"
echo "RUNNING ${hesin_script}" | tee -a "$logfile"
echo "" | tee -a "$logfile"

# create output file with samples having ICDs from exclusion criteria
PYTHONPATH="${upperdir}"/python "${hesin_script}" -p OA -r "$hesin_release" -o "${tempfiles[tmp_hesin_out]}" --icd9 "${tempfiles[tmp_icd9]}" --icd10 "${tempfiles[tmp_icd10]}" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)

n=$(cat "${tempfiles[tmp_hesin_out]}"| wc -l)
echo "" | tee -a "$logfile"
echo "INFO: $n samples meet exclusion criteria" | tee -a "$logfile"
echo "" | tee -a "$logfile"

# final output: from the list of samples meeting inclusion criteria, remove samples meeting exclusion criteria
join -1 1 -2 1 -a 1 -t$'\t' -e NULL -o 1.1,2.1 <(sort "${tempfiles[tmp_ukbb_out]}") <(sort "${tempfiles[tmp_hesin_out]}")|grep NULL|cut -f 1 >"$outfile"
n=$(cat "$outfile"|wc -l)
echo "INFO: final set: $n samples" | tee -a "$logfile"
echo "" | tee -a "$logfile"

# delete temp files
for fn in "${tempfiles[@]}";do
    if [[ -f "$fn" ]];then
	rm "$fn"
    fi
done

date "+%d-%b-%Y:%H-%M-%S"|tee -a "$logfile"

exit 0
