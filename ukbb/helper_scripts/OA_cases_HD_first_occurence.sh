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
    echo "                       -p | --keep <optional: keep temporary files>"
    echo "                       -h | --help"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

OPTS=$(getopt -o hpe:c:r:o:k: -l help,keep,hesin:,release:,config:,output:,key: -n 'OA_select_HD_first_occurence' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

declare -A keys=( [FingerOA]=1 [HandOA]=1 [HipKneeOA]=1 [HipOA]=1 [KneeOA]=1 [OA]=1 [SpineOA]=1 [ThumbOA]=1 )

declare -a tempfiles
config=""
hesin_release=""
outprefix=""
key=""
outfile=""
keep="NO"
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
    -p|--keep ) keep="YES"; shift ;;
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

tmpdir=$(mktemp -d -p "$out_dir" cases_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temp dir in $out_dir" | tee -a "$logfile"
    exit 1
else
    echo "INFO: temp dir: $tmpdir"    
fi

#----------------------------------------------------------------------------------------------------------------

date "+%d-%b-%Y:%H-%M-%S"|tee -a "$logfile"
echo ""|tee -a "$logfile"
echo "Current dir: ${PWD}"|tee -a "$logfile"
echo "Command line: $scriptname ${args[@]}"|tee -a "$logfile"
echo ""|tee -a "$logfile"
echo "INFO: HESIN release: $hesin_release"|tee -a "$logfile"
echo "INFO: key: $key"|tee -a "$logfile"
echo "INFO: temp dir: $tmpdir"|tee -a "$logfile"
echo "INFO: keep temp files: $keep"|tee -a "$logfile"
echo "INFO: output prefix: $outprefix"|tee -a "$logfile"
echo "INFO: output file: $outfile"|tee -a "$logfile"
echo ""|tee -a "$logfile"

#----------------------------------------------------------------------------------------------------------------

# inclusion ICD9 codes

while read code;do
    "$date_script" --icd9 "$code" -p "OA" -r "$hesin_release" -o "$tmpdir"/inclusion_icd9_"$code" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
    if [[ -f "$tmpdir"/inclusion_icd9_"$code" ]];then
	tempfiles+=("$tmpdir"/inclusion_icd9_"$code")
    else
	echo "DEBUG: no such file" "$tmpdir"/inclusion_icd9_"$code"
    fi    
done < <(grep -v "^#" "$icd_inclusion_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd9"{print $3;}')

# grep -v "^#" "$icd_inclusion_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd9"{print $3;}'|while read code;do
#     "$date_script" --icd9 "$code" -p "OA" -r "$hesin_release" -o "$tmpdir"/inclusion_icd9_"$code" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
#     if [[ -f "$tmpdir"/inclusion_icd9_"$code" ]];then
# 	tempfiles+=("$tmpdir"/inclusion_icd9_"$code")
#     else
# 	echo "DEBUG: no such file" "$tmpdir"/inclusion_icd9_"$code"
#     fi
# done

# inclusion ICD10 codes
while read code;do
    "$date_script" --icd10 "$code" -p "OA" -r "$hesin_release" -o "$tmpdir"/inclusion_icd10_"$code" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
    if [[ -f "$tmpdir"/inclusion_icd10_"$code" ]];then
	tempfiles+=("$tmpdir"/inclusion_icd10_"$code")
    else
	echo "DEBUG: no such file" "$tmpdir"/inclusion_icd10_"$code"
    fi
done < <(grep -v "^#" "$icd_inclusion_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd10"{print $3;}')

# grep -v "^#" "$icd_inclusion_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd10"{print $3;}'|while read code;do
#     "$date_script" --icd10 "$code" -p "OA" -r "$hesin_release" -o "$tmpdir"/inclusion_icd10_"$code" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
#     if [[ -f "$tmpdir"/inclusion_icd10_"$code" ]];then
# 	tempfiles+=("$tmpdir"/inclusion_icd10_"$code")
#     else
# 	echo "DEBUG: no such file" "$tmpdir"/inclusion_icd10_"$code"
#     fi
# done

echo "DEBUG: # tempfiles: ${#tempfiles[@]}"
merged_incl_fname="$tmpdir"/inclusion_merged
: > "$merged_incl_fname"
for f in "${tempfiles[@]}";do
    join_two_files "$merged_incl_fname" 1 "$f" 1 "$merged_incl_fname"
done	

tail -n +2 "$merged_incl_fname" | perl -lne 'BEGIN{$,="\t";sub compare{return -1 if $b!~/\d{2}\/\d{2}\/\d{4}/;return 1 if $a!~/\d{2}\/\d{2}\/\d{4}/;@m1=$a=~m/(\d{2})\/(\d{2})\/(\d{4})/;@m2=$b=~m/(\d{2})\/(\d{2})\/(\d{4})/;return -1 if ($m1[2]<$m2[2]);return 1 if ($m1[2]>$m2[2]);return -1 if ($m1[1]<$m2[1]);return 1 if ($m1[1]>$m2[1]);return -1 if ($m1[0]<$m2[0]);return 1 if ($m1[0]>$m2[0]);return 0;}}{@a=split(/\t/);$id=shift(@a);@b=sort compare @a;print $id,$b[0];}' > "$tmpdir"/inclusion_final

#------------------------------------------------------------------------------------------------------------------

tempfiles=()
# exclusion ICD9 codes
while read code;do
    "$status_script" --icd9 "$code" -p "OA" -r "$hesin_release" -o "$tmpdir"/exclusion_icd9_"$code" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
    if [[ -f "$tmpdir"/exclusion_icd9_"$code" ]];then    
	tempfiles+=("$tmpdir"/exclusion_icd9_"$code")
    fi
done < <(grep ^icd9 "$icd_exclusion_file"| cut -f 2)

# grep ^icd9 "$icd_exclusion_file"| cut -f 2| while read code;do
#     "$status_script" --icd9 "$code" -p "OA" -r "$hesin_release" -o "$tmpdir"/exclusion_icd9_"$code" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
#     if [[ -f "$tmpdir"/exclusion_icd9_"$code" ]];then    
# 	tempfiles+=("$tmpdir"/exclusion_icd9_"$code")
#     fi
# done

# exclusion ICD10 codes
while read code;do
    "$status_script" --icd10 "$code" -p "OA" -r "$hesin_release" -o "$tmpdir"/exclusion_icd10_"$code" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
    if [[ -f "$tmpdir"/exclusion_icd10_"$code" ]];then    
	tempfiles+=("$tmpdir"/exclusion_icd10_"$code")
    fi
done < <(grep ^icd10 "$icd_exclusion_file"| cut -f 2)

# grep ^icd10 "$icd_exclusion_file"| cut -f 2| while read code;do
#     "$status_script" --icd10 "$code" -p "OA" -r "$hesin_release" -o "$tmpdir"/exclusion_icd10_"$code" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
#     if [[ -f "$tmpdir"/exclusion_icd10_"$code" ]];then    
# 	tempfiles+=("$tmpdir"/exclusion_icd10_"$code")
#     fi
# done

echo "DEBUG: # tempfiles: ${#tempfiles[@]}"
merged_excl_fname="$tmpdir"/exclusion_merged
: > "$merged_excl_fname"
for f in "${tempfiles[@]}";do
    join_two_files "$merged_excl_fname" 1 "$f" 1 "$merged_excl_fname"
done	

tail -n +2 "$merged_excl_fname" | perl -lne 'BEGIN{$,="\t";}{@a=split(/\t/);$id=shift(@a);$s=0;foreach (@a){$s=1 if $_ == 1;}print $id,$s;}' > "$tmpdir"/exclusion_final

#------------------------------------------------------------------------------------------------------------------

join -t $'\t' -1 1 -2 1 -a 1 -a 2 -e "NA" -o 1.1,2.1,1.2,2.2 <(sort -k1,1 "$tmpdir"/inclusion_final) <(sort -k1,1 "$tmpdir"/exclusion_final) | awk 'BEGIN{FS=OFS="\t";}{if ($2=="NA"){$2=$1;}if ($4=="0"){print $2,$3;}}' > "$outfile"

#------------------------------------------------------------------------------------------------------------------

if [[ "$keep" == "NO" ]];then
    echo "INFO: deleting temp dir $tmpdir"
    rm -rf "$tmpdir"
fi

date "+%d-%b-%Y:%H-%M-%S"|tee -a "$logfile"

exit 0



