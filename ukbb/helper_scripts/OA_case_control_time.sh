#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
upperdir=$(dirname $scriptdir)
source "${upperdir}/functions.sh"
date_script="${upperdir}/first_diagnosis_date.py"
main_select_script="${upperdir}/ukbb_select.py"

function usage () {
    echo ""
    echo "Script for selecting OA case/control individuals for a specific time point"
    echo ""
    echo "Usage: OA_select_HD.sh -e | --hesin <release of the HESIN dataset>"
    echo "                       -r | --release <release of the MAIN dataset>"
    echo "                       -i | --instance <time point: 0,1,2,3>"
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

OPTS=$(getopt -o hpc:r:o:e:i: -l help,keep,config:,release:,output:,hesin:,instance: -n 'OA_case_control_time' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

declare -A oa_keys=( [FingerOA]=1 [HandOA]=1 [HipKneeOA]=1 [HipOA]=1 [KneeOA]=1 [OA]=1 [SpineOA]=1 [ThumbOA]=1 )
declare -a tmp_ar

OA_case_value="1465"

# declare -a tempfiles
config=""
hesin_release=""
main_release=""
outprefix=""
instance=""
outfile=""
keep="NO"
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
    -p|--keep ) keep="YES"; shift ;;
    -e|--hesin ) hesin_release=$2; shift 2 ;;
    -r|--release ) main_release=$2; shift 2 ;;
    -i|--instance ) instance=$2; shift 2 ;;
    -c|--config ) config=$2; shift 2 ;;
    -o|--output ) outprefix=$2; shift 2 ;;
    --) shift ; break ;;
    * ) echo "Unexpected option: $1" ; usage ;;
  esac
done

exitIfEmpty "$hesin_release" "ERROR: HESIN release not specified"
exitIfEmpty "$main_release" "ERROR: MAIN release not specified"
exitIfEmpty "$instance" "ERROR: instance not specified"
if [[ "$instance" -lt 0 || "$instance" -gt 3 ]];then
    echo "ERROR: instance (-i) values may only be 0, 1, 2, 3"
    exit 1
fi

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

op_inclusion_file=""
readValue "$config" HD_OA_OPCS4 op_inclusion_file
exitIfNotFile "$op_inclusion_file" "ERROR: HD_OA_OPCS4 ($op_inclusion_file) does not exist"

logfile="$outprefix".log
outfile_case="$outprefix".case.txt
outfile_control="$outprefix".control.txt
exitIfExists "$outfile_case" "ERROR: output file $outfile_case already exists"
exitIfExists "$outfile_control" "ERROR: output file $outfile_control already exists"
: > "$logfile"

out_dir=$(dirname "$outfile_case")

#----------------------------------------------------------------------------------------------------------------

tmpdir=$(mktemp -d -p "$out_dir" case_control_XXXXXXXX)
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
echo "INFO: MAIN release: $main_release"|tee -a "$logfile"
echo "INFO: HESIN release: $hesin_release"|tee -a "$logfile"
echo "INFO: temp dir: $tmpdir"|tee -a "$logfile"
echo "INFO: keep temp files: $keep"|tee -a "$logfile"
echo "INFO: output prefix: $outprefix"|tee -a "$logfile"
echo "INFO: output file: $outfile"|tee -a "$logfile"
echo ""|tee -a "$logfile"

#----------------------------------------------------------------------------------------------------------------
# 01
# SR cases, MAIN dataset
# defining case as someone having 1465 in field 20002 at the specified instance or earlier

tmpout="$tmpdir"/01.temp_main_20002.txt
"$main_select_script" -p "OA" -r "$main_release" -f "20002" -o "$tmpout"  > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
tmp_ar=(1) # ID and columns for instances <= given instance
while read c;do
    tmp_ar+=($c)
done < <(head -n 1 "$tmpout" | perl -slne '@a=split(/\t/);for ($i=0;$i<scalar(@a);$i++){$z=$a[$i];if ($z=~m/^f\.\d+\.(\d+)\.\d+/){print $i+1 if ($1<=$x);}}' -- -x="$instance")
oa_main_out="$tmpdir"/01.oa_main.txt
cut -f $(join_by "," "${tmp_ar[@]}") "$tmpout" | tail -n +2 | perl -slne '@a=split(/\t/);for ($i=1;$i<scalar(@a);$i++){if ($a[$i] eq $v){print $a[0];next;}}' -- -v="$OA_case_value" > "$oa_main_out"

#----------------------------------------------------------------------------------------------------------------
# 02
# VISIT DATES

visit_date="$tmpdir"/02.visits.txt
"$main_select_script" -p "OA" -r "$main_release" -f "53" -o "$visit_date"  > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
instcol=$(head -n 1 "$visit_date" | perl -slne '@a=split(/\t/);for ($i=0;$i<scalar(@a);$i++){if ($a[$i]=~m/^f\.\d+\.(\d)\.\d+/){print $i+1 if ($1==$x);}}' -- -x="$instance")
cut -f 1,"$instcol" "$visit_date" | tail -n +2 | grep -v "NA" | perl -lne '@a=split(/\t/);if ($a[1]=~m/(\d{4})-(\d{2})-(\d{2})/){$a[1]=$3."/".$2."/".$1;} print join("\t",@a);' | sponge "$visit_date" # ID and given instance columns, no header, skip NAs

#----------------------------------------------------------------------------------------------------------------
# 03
# inclusion ICD codes

declare -a inclusion_final_files
for key in "${!oa_keys[@]}";do
    tmp_incl_icd9="$tmpdir"/03.inclusion_"${key}"_icd9
    tmp_incl_icd10="$tmpdir"/03.inclusion_"${key}"_icd10
    tmp_incl_merged="$tmpdir"/03.inclusion_merged_"${key}"
    tmp_incl_final="$tmpdir"/03.inclusion_"${key}"_final
    
    "$date_script" --icd9 <(grep -v "^#" "$icd_inclusion_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd9"{print $3;}') -p "OA" -r "$hesin_release" -o "$tmp_incl_icd9" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
    "$date_script" --icd10 <(grep -v "^#" "$icd_inclusion_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd10"{print $3;}') -p "OA" -r "$hesin_release" -o "$tmp_incl_icd10" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
    
    if [[ -f "$tmp_incl_icd9" && -f "$tmp_incl_icd10" ]];then
	join_two_files "$tmp_incl_icd9" 1 "$tmp_incl_icd10" 1 "$tmp_incl_merged" # output with header
	tail -n +2 "$tmp_incl_merged" | perl -lne 'BEGIN{$,="\t";sub compare{return -1 if $b!~/\d{2}\/\d{2}\/\d{4}/;return 1 if $a!~/\d{2}\/\d{2}\/\d{4}/;@m1=$a=~m/(\d{2})\/(\d{2})\/(\d{4})/;@m2=$b=~m/(\d{2})\/(\d{2})\/(\d{4})/;return -1 if ($m1[2]<$m2[2]);return 1 if ($m1[2]>$m2[2]);return -1 if ($m1[1]<$m2[1]);return 1 if ($m1[1]>$m2[1]);return -1 if ($m1[0]<$m2[0]);return 1 if ($m1[0]>$m2[0]);return 0;}}{@a=split(/\t/);$id=shift(@a);@b=sort compare @a;print $id,$b[0] if $b[0]=~/\d{2}\/\d{2}\/\d{4}/;}' > "$tmp_incl_final" # no header
	inclusion_final_files+=("$tmp_incl_final")
    fi
done

inclusion_final="$tmpdir"/03.inclusion_allkeys_final
# join with visit dates, compare two dates, output sample ID if the first date <= second date
for f in "${inclusion_final_files[@]}";do
    join -1 1 -2 1 -t$'\t' -o 1.1,1.2,2.2 <(sort -k1,1 "$f") <(sort -k1,1 "$visit_date") | perl -lne '@a=split(/\t/);$id=shift(@a);@m1=$a[0]=~m/(\d{2})\/(\d{2})\/(\d{4})/;@m2=$a[1]=~m/(\d{2})\/(\d{2})\/(\d{4})/;next if ($m1[2]>$m2[2]);if ($m1[2]<$m2[2]){print $id;next;}next if ($m1[1]>$m2[1]);if ($m1[1]<$m2[1]){print $id;next;}next if ($m1[2]>$m2[2]);if ($m1[2]<=$m2[2]){print $id;}'
done | sort | uniq > "$inclusion_final"

#------------------------------------------------------------------------------------------------------------------
# 04
# union of all cases

all_cases="$tmpdir"/04.all_cases
cat "$inclusion_final" "$oa_main_out" | sort | uniq > "$all_cases"

#------------------------------------------------------------------------------------------------------------------
# 05
# exclusion ICD codes

temp_excl_icd9="$tmpdir"/05.exclusion_icd9
temp_excl_icd10="$tmpdir"/05.exclusion_icd10

"$date_script" --icd9 <(grep ^icd9 "$icd_exclusion_file" | cut -f 2) -p "OA" -r "$hesin_release" -o "$temp_excl_icd9" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)
"$date_script" --icd10 <(grep ^icd10 "$icd_exclusion_file" | cut -f 2) -p "OA" -r "$hesin_release" -o "$temp_excl_icd10" > >(tee -a "$logfile") 2> >(tee -a "$logfile" >&2)

join_two_files "$tmpdir"/exclusion_icd10 1 "$tmpdir"/exclusion_icd9 1 "$tmpdir"/exclusion_merged
tail -n +2 "$tmpdir"/exclusion_merged | perl -lne 'BEGIN{$,="\t";}{@a=split(/\t/);$id=shift(@a);$s=0;foreach (@a){$s=1 if $_ == 1;}print $id,$s;}' > "$tmpdir"/exclusion_final

#------------------------------------------------------------------------------------------------------------------

join -t $'\t' -1 1 -2 1 -a 1 -a 2 -e "NA" -o 1.1,2.1,1.2,2.2 <(sort -k1,1 "$tmpdir"/inclusion_final) <(sort -k1,1 "$tmpdir"/exclusion_final) | awk 'BEGIN{FS=OFS="\t";}{if ($2=="NA"){$2=$1;}if ($4=="0" && $3!="NA"){print $2,$3;}}' > "$outfile"

#------------------------------------------------------------------------------------------------------------------

if [[ "$keep" == "NO" ]];then
    echo "INFO: deleting temp dir $tmpdir"
    rm -rf "$tmpdir"
fi

date "+%d-%b-%Y:%H-%M-%S" | tee -a "$logfile"

exit 0



