#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
upperdir=$(dirname $scriptdir)
source "${upperdir}/functions.sh"
hesin_script="${upperdir}/hesin_select.py"
fullsname=$(realpath $0)

function usage () {
    echo ""
    echo "Script for selecting hospital diagnosed OA knee/hip replacement cases from a HESIN release"
    echo ""
    echo "Usage: OA_select_HD_replacement.sh -e | --hesin <release of the HESIN dataset>"
    echo "                                   -k | --key <one of: THR, TKR, TJR>"
    echo "                                   -o | --output <output prefix>"
    echo "                                   -f | --force <overwrite existing output>"
    echo "                                   -c | --config <optional: config file; default: config.txt in script directory>"
    echo "                                   -h | --help"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

OPTS=$(getopt -o hfc:e:r:o:k: -l help,force,hesin:,release:,config:,output:,key: -n 'OA_select_HD_replacement' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

declare -A keys=( [THR]=HipOA [TKR]=KneeOA [TJR]=1 )
declare -A tempfiles

config=""
hesin_release=""
outprefix=""
outfile=""
key=""
force="NO"
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
    -e|--hesin) hesin_release=$2; shift 2 ;;
    -k|--key ) key=$2; shift 2 ;;
    -c|--config ) config=$2; shift 2 ;;
    -o|--output ) outprefix=$2; shift 2 ;;
    -f|--force ) force="YES"; shift ;;
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

opcs4_inclusion_file=""
readValue "$config" HD_OA_OPCS4 opcs4_inclusion_file
exitIfNotFile "$opcs4_inclusion_file" "ERROR: HD_OA_OPCS4 ($opcs4_inclusion_file) does not exist"

icd_exclusion_file=""
readValue "$config" OA_CASE_EXCLUSION icd_exclusion_file
exitIfNotFile "$icd_exclusion_file" "ERROR: OA_CASE_EXCLUSION ($icd_exclusion_file) does not exist"

logfile="$outprefix".log
outfile="$outprefix".txt
if [[ "$force" == "NO" ]];then
    exitIfExists "$outfile" "ERROR: output file $outfile already exists"
fi

: > "$logfile"
out_dir=$(dirname "$outfile")

# for hesin_select, free CPUs
ncpus=$(getFreeCPUs 0.5)

#----------------------------------------------------------------------------------------------------------------

date "+%d-%b-%Y:%H-%M-%S"|tee -a "$logfile"
echo ""|tee -a "$logfile"
echo "Current dir: ${PWD}"|tee -a "$logfile"
echo "Command line: $scriptname ${args[@]}"|tee -a "$logfile"
echo ""|tee -a "$logfile"

echo "INFO: HESIN release: $hesin_release"|tee -a "$logfile"
echo "INFO: key: $key"|tee -a "$logfile"
echo "INFO: force: $force"|tee -a "$logfile"
echo "INFO: output prefix: $outprefix"|tee -a "$logfile"
echo "INFO: output file: $outfile"|tee -a "$logfile"
echo ""|tee -a "$logfile"

#----------------------------------------------------------------------------------------------------------------
if [[ "$key" == "TJR" ]];then # call itself with TKR and THR keys and take the union
    tempfiles["tmp_TKR_output"]="tmp_TKR_XXXXXXX.txt"
    tempfiles["tmp_THR_output"]="tmp_THR_XXXXXXX.txt"
    if createTempFiles "$out_dir" tempfiles;then
	echo "INFO: created temporary files"|tee -a "$logfile"
    else
	echo "ERROR: failed to create temporary files"|tee -a "$logfile"
	exit 1
    fi

    # TKR
    echo "RUNNING TKR" |tee -a "$logfile"
    prfx="${tempfiles[tmp_TKR_output]}"
    prfx=${prfx/%.txt}
    "$fullsname" -f -k TKR -e "$hesin_release" -o "$prfx" -c "${config}"
    cat "$prfx".log >> "$logfile"
    rm "$prfx".log
    echo "" |tee -a "$logfile"

    # THR
    echo "RUNNING THR" |tee -a "$logfile"
    prfx="${tempfiles[tmp_THR_output]}"
    prfx=${prfx/%.txt}
    "$fullsname" -f -k THR -e "$hesin_release" -o "$prfx" -c "${config}"
    cat "$prfx".log >> "$logfile"
    rm "$prfx".log
    echo "" |tee -a "$logfile"

    # take the union
    cat "${tempfiles[tmp_TKR_output]}" "${tempfiles[tmp_THR_output]}" | sort | uniq > "${outfile}"
    echo "INFO: final output list:" $(cat "${outfile}"|wc -l) "IDs"  | tee -a "$logfile"
    echo "" |tee -a "$logfile"
    
    # delete temp files
    for fn in "${tempfiles[@]}";do
	if [[ -f "$fn" ]];then
	    rm "$fn"
	fi
    done
    
    date "+%d-%b-%Y:%H-%M-%S"|tee -a "$logfile"
    exit 0
fi
#----------------------------------------------------------------------------------------------------------------

# create temp files
tempfiles["tmp_icd9_incl"]="tmp_XXXXXXX"
tempfiles["tmp_icd10_incl"]="tmp_XXXXXXX"
tempfiles["tmp_inclusion_icd_out"]="tmp_XXXXXXX"
tempfiles["tmp_opcs4_incl"]="tmp_XXXXXXX"
tempfiles["tmp_inclusion_opcs4_out"]="tmp_XXXXXXX"
tempfiles["tmp_inclusion_out"]="tmp_XXXXXXX"
tempfiles["tmp_icd9_excl"]="tmp_XXXXXXX"
tempfiles["tmp_icd10_excl"]="tmp_XXXXXXX"
tempfiles["tmp_exclusion_out"]="tmp_XXXXXXX"

if createTempFiles "$out_dir" tempfiles;then
    echo "INFO: created temporary files"|tee -a "$logfile"
else
    echo "ERROR: failed to create temporary files"|tee -a "$logfile"
    exit 1
fi

#----------------------------------------------------------------------------------------------------------------

# CASE INCLUSION
# inclusion ICD codes
grep -v "^#" "$icd_inclusion_file" | awk -v k="${keys[$key]}" 'BEGIN{FS="\t";}$1==k && $2=="icd9"{print $3;}'|sort|uniq > "${tempfiles[tmp_icd9_incl]}"
grep -v "^#" "$icd_inclusion_file" | awk -v k="${keys[$key]}" 'BEGIN{FS="\t";}$1==k && $2=="icd10"{print $3;}'|sort|uniq > "${tempfiles[tmp_icd10_incl]}"

echo ""  | tee -a "$logfile"
PYTHONPATH="${upperdir}"/python "$hesin_script" -p "OA" -r "$hesin_release" --icd9 "${tempfiles[tmp_icd9_incl]}" --icd10 "${tempfiles[tmp_icd10_incl]}" -o "${tempfiles[tmp_inclusion_icd_out]}" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
echo ""  | tee -a "$logfile"

# inclusion OPCS4 codes
if [[ "$key" == "TJR" ]];then # take all codes for both TKR and THR
    grep -v "^#" "$opcs4_inclusion_file" | cut -f 2| sort|uniq > "${tempfiles[tmp_opcs4_incl]}"
else # select based on key
    grep -v "^#" "$opcs4_inclusion_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k{print $2;}'| sort|uniq > "${tempfiles[tmp_opcs4_incl]}"
fi
echo ""  | tee -a "$logfile"
PYTHONPATH="${upperdir}"/python "$hesin_script" -p "OA" -r "$hesin_release" -n $ncpus --oper4 "${tempfiles[tmp_opcs4_incl]}" -o "${tempfiles[tmp_inclusion_opcs4_out]}" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
echo ""  | tee -a "$logfile"

# final inclusion set: taking intersection of ICD inclusion and OPCS4 inclusion lists
cat "${tempfiles[tmp_inclusion_icd_out]}" "${tempfiles[tmp_inclusion_opcs4_out]}"|sort|uniq -d > "${tempfiles[tmp_inclusion_out]}"
echo "INFO: case inclusion:" $(cat "${tempfiles[tmp_inclusion_out]}"| wc -l) "IDs"  | tee -a "$logfile"
echo ""  | tee -a "$logfile"

# CASE EXCLUSION
# exclusion ICD codes
grep ^icd9 "$icd_exclusion_file"| cut -f 2 > "${tempfiles[tmp_icd9_excl]}"
grep ^icd10 "$icd_exclusion_file"| cut -f 2 > "${tempfiles[tmp_icd10_excl]}"
PYTHONPATH="${upperdir}"/python "$hesin_script" -p "OA" -r "$hesin_release" --icd9 "${tempfiles[tmp_icd9_excl]}" --icd10 "${tempfiles[tmp_icd10_excl]}" -o "${tempfiles[tmp_exclusion_out]}" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
echo ""  | tee -a "$logfile"

# final output: only output samples in inclusion that are not in exclusion
join -1 1 -2 1 -a 1 -t$'\t' -e NULL -o 1.1,2.1 <(sort "${tempfiles[tmp_inclusion_out]}") <(sort "${tempfiles[tmp_exclusion_out]}") | grep NULL | cut -f 1  > "${outfile}"
echo "INFO: final output list:" $(cat "${outfile}"| wc -l) "IDs"  | tee -a "$logfile"
echo ""  | tee -a "$logfile"

# delete temp files
for fn in "${tempfiles[@]}";do
    if [[ -f "$fn" ]];then
	rm "$fn"
    fi
done

date "+%d-%b-%Y:%H-%M-%S"|tee -a "$logfile"
exit 0
