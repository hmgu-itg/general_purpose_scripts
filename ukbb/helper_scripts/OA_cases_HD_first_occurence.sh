#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
upperdir=$(dirname $scriptdir)
source "${upperdir}/functions.sh"
date_script="${upperdir}/first_diagnosis_date.py"
status_script="${upperdir}/get_icd_status.py"

function usage () {
    echo ""
    echo "Script for selecting first date of hospital diagnosed OA cases from a HESIN release"
    echo ""
    echo "Usage: OA_select_HD.sh -e | --hesin <release of the HESIN dataset>"
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
hesin_release=""
outprefix=""
key=""
outfile=""
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
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
exitIfEmpty "$hesin_release" "ERROR: HESIN release not specified"
exitIfEmpty "$outprefix" "ERROR: output prefix not specified"
if [[ -z "$config" ]];then
    config="${upperdir}"/config.txt
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

declare -a tempfiles
# tempfiles["tmp_icd9_incl"]="tmp_XXXXXXX"
# tempfiles["tmp_icd10_incl"]="tmp_XXXXXXX"
# tempfiles["tmp_inclusion_out"]="tmp_XXXXXXX"
# tempfiles["tmp_icd9_excl"]="tmp_XXXXXXX"
# tempfiles["tmp_icd10_excl"]="tmp_XXXXXXX"
# tempfiles["tmp_exclusion_out"]="tmp_XXXXXXX"

# if createTempFiles "$out_dir" tempfiles;then
#     echo "INFO: created temporary files"|tee -a "$logfile"
# else
#     echo "ERROR: failed to create temporary files"|tee -a "$logfile"
#     exit 1
# fi

tmpdir=$(mktemp -d -p "$out_dir" cases_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temp dir in $out_dir" | tee -a "$logfile"
    exit 1
fi

#----------------------------------------------------------------------------------------------------------------

date "+%d-%b-%Y:%H-%M-%S"|tee -a "$logfile"
echo ""|tee -a "$logfile"
echo "Current dir: ${PWD}"|tee -a "$logfile"
echo "Command line: $scriptname ${args[@]}"|tee -a "$logfile"
echo ""|tee -a "$logfile"
echo "INFO: HESIN release: $hesin_release"|tee -a "$logfile"
echo "INFO: key: $key"|tee -a "$logfile"
echo "INFO: output prefix: $outprefix"|tee -a "$logfile"
echo "INFO: output file: $outfile"|tee -a "$logfile"
echo ""|tee -a "$logfile"

#----------------------------------------------------------------------------------------------------------------

# inclusion ICD9 codes
grep -v "^#" "$icd_inclusion_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd9"{print $3;}'|while read code;do
    "$date_script" --icd9 "$code" -p "OA" -r "$hesin_release" -o "$tmpdir"/inclusion_icd9_"$code" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
    tempfiles+=("$tmpdir"/inclusion_icd9_"$code")
done

# inclusion ICD10 codes
grep -v "^#" "$icd_inclusion_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd10"{print $3;}'|while read code;do
    "$date_script" --icd10 "$code" -p "OA" -r "$hesin_release" -o "$tmpdir"/inclusion_icd10_"$code" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
    tempfiles+=("$tmpdir"/inclusion_icd10_"$code")
done

merged_incl_fname="$tmpdir"/inclusion_merged
: > "$merged_incl_fname"
for f in "${tempfiles[@]}";do
    join_two_files "$merged_incl_fname" 1 "$f" 1 "$merged_incl_fname"
done	

tail -n +2 "$merged_incl_fname" | perl -lne 'BEGIN{$,="\t";sub compare{return -1 if $b!~/\d{2}\/\d{2}\/\d{4}/;return 1 if $a!~/\d{2}\/\d{2}\/\d{4}/;@m1=$a=~m/(\d{2})\/(\d{2})\/(\d{4})/;@m2=$b=~m/(\d{2})\/(\d{2})\/(\d{4})/;return -1 if ($m1[2]<$m2[2]);return 1 if ($m1[2]>$m2[2]);return -1 if ($m1[1]<$m2[1]);return 1 if ($m1[1]>$m2[1]);return -1 if ($m1[0]<$m2[0]);return 1 if ($m1[0]>$m2[0]);return 0;}}{@a=split(/\t/);$id=shift(@a);@b=sort compare @a;print $id,$b[0];}' > "$tmpdir"/inclusion_final

#------------------------------------------------------------------------------------------------------------------

tempfiles=()
# exclusion ICD9 codes
grep ^icd9 "$icd_exclusion_file"| cut -f 2| while read code;do
    "$status_script" --icd9 "$code" -p "OA" -r "$hesin_release" -o "$tmpdir"/exclusion_icd9_"$code"
    tempfiles+=("$tmpdir"/exclusion_icd9_"$code")
done

# exclusion ICD10 codes
grep ^icd10 "$icd_exclusion_file"| cut -f 2| while read code;do
    "$status_script" --icd10 "$code" -p "OA" -r "$hesin_release" -o "$tmpdir"/exclusion_icd10_"$code"
    tempfiles+=("$tmpdir"/exclusion_icd10_"$code")
done

merged_excl_fname="$tmpdir"/exclusion_merged
: > "$merged_excl_fname"
for f in "${tempfiles[@]}";do
    join_two_files "$merged_excl_fname" 1 "$f" 1 "$merged_excl_fname"
done	

tail -n +2 "$merged_excl_fname" | perl -lne 'BEGIN{$,="\t";}{@a=split(/\t/);$id=shift(@a);$s=0;foreach (@a){$s=1 if $_ == 1;}print $id,$s;}' > "$tmpdir"/exclusion_final

#------------------------------------------------------------------------------------------------------------------




grep ^icd9 "$icd_exclusion_file"| cut -f 2 > "${tempfiles[tmp_icd9_excl]}"
grep ^icd10 "$icd_exclusion_file"| cut -f 2 > "${tempfiles[tmp_icd10_excl]}"
echo ""  | tee -a "$logfile"
PYTHONPATH="${upperdir}"/python "$hesin_script" -p "OA" -r "$hesin_release" --icd9 "${tempfiles[tmp_icd9_excl]}" --icd10 "${tempfiles[tmp_icd10_excl]}" -o "${tempfiles[tmp_exclusion_out]}" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)

# only output samples in inclusion that are not in exclusion
join -1 1 -2 1 -a 1 -t$'\t' -e NULL -o 1.1,2.1 <(sort "${tempfiles[tmp_inclusion_out]}") <(sort "${tempfiles[tmp_exclusion_out]}") | grep NULL | cut -f 1  > "${outfile}"
echo ""  | tee -a "$logfile"
echo "INFO: final output list:" $(cat "${outfile}"| wc -l) "IDs" | tee -a "$logfile"
echo ""  | tee -a "$logfile"

# delete temp files
for fn in "${tempfiles[@]}";do
    if [[ -f "$fn" ]];then
	rm "$fn"
    fi
done

date "+%d-%b-%Y:%H-%M-%S"|tee -a "$logfile"

exit 0



