#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
upperdir=$(dirname $scriptdir)
source "${upperdir}/functions.sh"
filter_script="${upperdir}/filter_by_visit.pl"
date_script="${upperdir}/first_diagnosis_date.py"
date_op_script="${upperdir}/first_operation_date.py"
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
exitIfDir "$outprefix" "ERROR: output prefix is a directory"
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

icd_control_exclusion_file=""
readValue "$config" OA_CONTROL_EXCLUSION_ICD icd_control_exclusion_file
exitIfNotFile "$icd_control_exclusion_file" "ERROR: OA_CONTROL_EXCLUSION_ICD ($icd_control_exclusion_file) does not exist"

op_control_exclusion_file=""
readValue "$config" OA_CONTROL_EXCLUSION_OPCS4 op_control_exclusion_file
exitIfNotFile "$op_control_exclusion_file" "ERROR: OA_CONTROL_EXCLUSION_OPCS4 ($op_control_exclusion_file) does not exist"

op_inclusion_file=""
readValue "$config" HD_OA_OPCS4 op_inclusion_file
exitIfNotFile "$op_inclusion_file" "ERROR: HD_OA_OPCS4 ($op_inclusion_file) does not exist"

logfile="$outprefix".log
outfile="$outprefix".txt
exitIfExists "$outfile" "ERROR: output file $outfile already exists"
: > "$logfile"

out_dir=$(dirname "$outfile")

#----------------------------------------------------------------------------------------------------------------

tmpdir=$(mktemp -d -p "$out_dir" case_control_time_XXXXXXXX)
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
echo "MAIN release: $main_release"|tee -a "$logfile"
echo "HESIN release: $hesin_release"|tee -a "$logfile"
echo "instance: $instance"|tee -a "$logfile"
echo "temp dir: $tmpdir"|tee -a "$logfile"
echo "keep temp files: $keep"|tee -a "$logfile"
echo "output prefix: $outprefix"|tee -a "$logfile"
echo "output file: $outfile"|tee -a "$logfile"
echo ""|tee -a "$logfile"

#----------------------------------------------------------------------------------------------------------------
# 01
# SR cases, MAIN dataset
# defining case as someone having 1465 in field 20002 at the specified instance or earlier

tmpout="$tmpdir"/01.1.main_20002.txt # with header
oa_main_out="$tmpdir"/01.2.oa_main_cases.txt # no header

echo $(date "+%d-%b-%Y:%H-%M-%S:") "01 selecting SR cases; saving results in $oa_main_out"

"$main_select_script" -p "OA" -r "$main_release" -f "20002" -o "$tmpout" >>"$logfile" 2>>"$logfile"
tmp_ar=(1) # ID and columns for instances <= given instance
while read c;do
    tmp_ar+=($c)
done < <(head -n 1 "$tmpout" | perl -slne '@a=split(/\t/);for ($i=0;$i<scalar(@a);$i++){$z=$a[$i];if ($z=~m/^f\.\d+\.(\d+)\.\d+/){print $i+1 if ($1<=$x);}}' -- -x="$instance")

cut -f $(join_by "," "${tmp_ar[@]}") "$tmpout" | tail -n +2 | perl -slne '@a=split(/\t/);for ($i=1;$i<scalar(@a);$i++){if ($a[$i] eq $v){print $a[0];next;}}' -- -v="$OA_case_value" > "$oa_main_out"

#----------------------------------------------------------------------------------------------------------------
# 02
# visit dates

visit_date="$tmpdir"/02.visits.txt # with header

echo $(date "+%d-%b-%Y:%H-%M-%S:") "02 selecting visit dates; saving results in $visit_date"

"$main_select_script" -p "OA" -r "$main_release" -f "53" -o "$visit_date" >>"$logfile" 2>>"$logfile"

# instcol=$(head -n 1 "$visit_date" | perl -slne '@a=split(/\t/);for ($i=0;$i<scalar(@a);$i++){if ($a[$i]=~m/^f\.\d+\.(\d)\.\d+/){print $i+1 if ($1==$x);}}' -- -x="$instance")
# echo "DEBUG: column for instance $instance: $instcol"
# cut -f 1,"$instcol" "$visit_date" | tail -n +2 | grep -v "NA" | perl -lne '@a=split(/\t/);if ($a[1]=~m/(\d{4})-(\d{2})-(\d{2})/){$a[1]=$3."/".$2."/".$1;} print join("\t",@a);' | sponge "$visit_date" # ID and given instance column, re-formatted date, no header, skip NAs

#----------------------------------------------------------------------------------------------------------------
# 03
# cases, inclusion ICD codes


inclusion_final="$tmpdir"/03.5.case_inclusion_allkeys_final # no header
echo $(date "+%d-%b-%Y:%H-%M-%S:") "03 selecting cases based on ICD codes; saving results in $inclusion_final"

declare -a inclusion_final_files
for key in "${!oa_keys[@]}";do
    tmp_incl_icd9="$tmpdir"/03.1.case_inclusion_"${key}"_icd9
    tmp_incl_icd10="$tmpdir"/03.2.case_inclusion_"${key}"_icd10
    tmp_incl_merged="$tmpdir"/03.3.case_inclusion_merged_"${key}"
    tmp_incl_min="$tmpdir"/03.4.case_inclusion_"${key}"_min
    
    "$date_script" --icd9 <(grep -v "^#" "$icd_inclusion_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd9"{print $3;}') -p "OA" -r "$hesin_release" -o "$tmp_incl_icd9" >>"$logfile" 2>>"$logfile"
    "$date_script" --icd10 <(grep -v "^#" "$icd_inclusion_file" | awk -v k=$key 'BEGIN{FS="\t";}$1==k && $2=="icd10"{print $3;}') -p "OA" -r "$hesin_release" -o "$tmp_incl_icd10" >>"$logfile" 2>>"$logfile"
    
    join_two_files "$tmp_incl_icd9" 1 "$tmp_incl_icd10" 1 "$tmp_incl_merged" # output with header
    
    if [[ -e "$tmp_incl_merged" ]];then
	tail -n +2 "$tmp_incl_merged" | perl -lne 'BEGIN{$,="\t";sub compare{return -1 if $b!~/\d{2}\/\d{2}\/\d{4}/;return 1 if $a!~/\d{2}\/\d{2}\/\d{4}/;@m1=$a=~m/(\d{2})\/(\d{2})\/(\d{4})/;@m2=$b=~m/(\d{2})\/(\d{2})\/(\d{4})/;return -1 if ($m1[2]<$m2[2]);return 1 if ($m1[2]>$m2[2]);return -1 if ($m1[1]<$m2[1]);return 1 if ($m1[1]>$m2[1]);return -1 if ($m1[0]<$m2[0]);return 1 if ($m1[0]>$m2[0]);return 0;}}{@a=split(/\t/);$id=shift(@a);@b=sort compare @a;print $id,$b[0];}' > "$tmp_incl_min" # no header
	inclusion_final_files+=("$tmp_incl_min")
    fi
done

if [[ "${#inclusion_final_files[@]}" -ne 0 ]];then
    for f in "${inclusion_final_files[@]}";do
	"$filter_script" "$f" "$visit_date" "$instance"
	# join -1 1 -2 1 -t$'\t' -o 1.1,1.2,2.2 <(sort -k1,1 "$f") <(sort -k1,1 "$visit_date") | perl -lne '@a=split(/\t/);$id=shift(@a);@m1=$a[0]=~m/(\d{2})\/(\d{2})\/(\d{4})/;@m2=$a[1]=~m/(\d{2})\/(\d{2})\/(\d{4})/;next if ($m1[2]>$m2[2]);if ($m1[2]<$m2[2]){print $id;next;}next if ($m1[1]>$m2[1]);if ($m1[1]<$m2[1]){print $id;next;}next if ($m1[0]>$m2[0]);if ($m1[0]<=$m2[0]){print $id;}'
    done | sort | uniq > "$inclusion_final"
else
    : > "$inclusion_final"
fi

#------------------------------------------------------------------------------------------------------------------
# 04
# union of all cases based on inclusion criteria

union_cases="$tmpdir"/04.union_cases # no header
echo $(date "+%d-%b-%Y:%H-%M-%S:") "04 combining SR cases and ICD based cases; saving results in $union_cases"
cat "$inclusion_final" "$oa_main_out" | sort | uniq > "$union_cases"

#------------------------------------------------------------------------------------------------------------------
# 05
# cases, exclusion ICD codes

temp_excl_icd9="$tmpdir"/05.1.exclusion_icd9
temp_excl_icd10="$tmpdir"/05.2.exclusion_icd10
temp_excl_merged="$tmpdir"/05.3.exclusion_merged
temp_excl_min="$tmpdir"/05.4.exclusion_min
temp_excl_final="$tmpdir"/05.5.exclusion_final

echo $(date "+%d-%b-%Y:%H-%M-%S:") "05 excluding cases based on ICD codes; saving results in $temp_excl_final"

"$date_script" --icd9 <(grep ^icd9 "$icd_exclusion_file" | cut -f 2) -p "OA" -r "$hesin_release" -o "$temp_excl_icd9" >>"$logfile" 2>>"$logfile"
"$date_script" --icd10 <(grep ^icd10 "$icd_exclusion_file" | cut -f 2) -p "OA" -r "$hesin_release" -o "$temp_excl_icd10" >>"$logfile" 2>>"$logfile"

join_two_files "$temp_excl_icd9" 1 "$temp_excl_icd10" 1 "$temp_excl_merged" # output with header

if [[ -e "$temp_excl_merged" ]];then
    tail -n +2 "$temp_excl_merged" | perl -lne 'BEGIN{$,="\t";sub compare{return -1 if $b!~/\d{2}\/\d{2}\/\d{4}/;return 1 if $a!~/\d{2}\/\d{2}\/\d{4}/;@m1=$a=~m/(\d{2})\/(\d{2})\/(\d{4})/;@m2=$b=~m/(\d{2})\/(\d{2})\/(\d{4})/;return -1 if ($m1[2]<$m2[2]);return 1 if ($m1[2]>$m2[2]);return -1 if ($m1[1]<$m2[1]);return 1 if ($m1[1]>$m2[1]);return -1 if ($m1[0]<$m2[0]);return 1 if ($m1[0]>$m2[0]);return 0;}}{@a=split(/\t/);$id=shift(@a);@b=sort compare @a;print $id,$b[0];}' > "$temp_excl_min" # no header
    "$filter_script" "$temp_excl_min" "$visit_date" "$instance" > "$temp_excl_final"
    # join with visit dates, compare two dates, output sample ID if the first date <= second date
    # join -1 1 -2 1 -t$'\t' -o 1.1,1.2,2.2 <(sort -k1,1 "$temp_excl_min") <(sort -k1,1 "$visit_date") | perl -lne '@a=split(/\t/);$id=shift(@a);@m1=$a[0]=~m/(\d{2})\/(\d{2})\/(\d{4})/;@m2=$a[1]=~m/(\d{2})\/(\d{2})\/(\d{4})/;next if ($m1[2]>$m2[2]);if ($m1[2]<$m2[2]){print $id;next;}next if ($m1[1]>$m2[1]);if ($m1[1]<$m2[1]){print $id;next;}next if ($m1[0]>$m2[0]);if ($m1[0]<=$m2[0]){print $id;}' > "$temp_excl_final"
fi

#------------------------------------------------------------------------------------------------------------------
# 06
# final cases

final_cases="$tmpdir"/06.final_cases # no header
echo $(date "+%d-%b-%Y:%H-%M-%S:") "06 creating final set of cases; saving results in $final_cases"

if [[ -e "$temp_excl_final" ]];then
    cat "$union_cases" | perl -snle 'BEGIN{%h=();if (length($f)!=0){open(fh,"<",$f);while(<fh>){chomp;$h{$_}=1;}close(fh);}}{@a=split(/\t/);if (!defined($h{$a[0]})){print $_;}}' -- -f="$temp_excl_final" > "$final_cases"
    # join -t $'\t' -1 1 -2 1 -a 1 -a 2 -e "NA" -o 1.1,2.1 <(sort "$union_cases") <(sort "$temp_excl_final") | awk 'BEGIN{FS=OFS="\t";}{if ($2=="NA"){print $1;}}' > "$final_cases"
else
    cp "$union_cases" "$final_cases"
fi

#------------------------------------------------------------------------------------------------------------------
# 07
# total samples

total_samples="$tmpdir"/07.total_samples # with header
echo $(date "+%d-%b-%Y:%H-%M-%S:") "07 extracting all sample IDs; saving results in $total_samples"

"$main_select_script" -p "OA" -r "$main_release" -o "$total_samples" >>"$logfile" 2>>"$logfile"

#------------------------------------------------------------------------------------------------------------------
# 08
# control exclusion ICD codes

ctl_excl_icd9="$tmpdir"/08.1.ctl_exclusion_icd9 # with header
ctl_excl_icd10="$tmpdir"/08.2.ctl_exclusion_icd10 # with header
ctl_excl_merged="$tmpdir"/08.3.ctl_exclusion_merged # with header
ctl_excl_min="$tmpdir"/08.4.ctl_exclusion_min # no header
ctl_excl_final_icd="$tmpdir"/08.5.ctl_exclusion_final_icd # no header

echo $(date "+%d-%b-%Y:%H-%M-%S:") "08 excluding samples from controls based on ICD codes; saving results in $ctl_excl_final_icd"

"$date_script" --icd9 <(grep ^icd9 "$icd_control_exclusion_file" | cut -f 2) -p "OA" -r "$hesin_release" -o "$ctl_excl_icd9" >>"$logfile" 2>>"$logfile"
"$date_script" --icd10 <(grep ^icd10 "$icd_control_exclusion_file" | cut -f 2) -p "OA" -r "$hesin_release" -o "$ctl_excl_icd10" >>"$logfile" 2>>"$logfile"
    
join_two_files "$ctl_excl_icd9" 1 "$ctl_excl_icd10" 1 "$ctl_excl_merged" # output with header

if [[ -e "$ctl_excl_merged" ]];then
    tail -n +2 "$ctl_excl_merged" | perl -lne 'BEGIN{$,="\t";sub compare{return -1 if $b!~/\d{2}\/\d{2}\/\d{4}/;return 1 if $a!~/\d{2}\/\d{2}\/\d{4}/;@m1=$a=~m/(\d{2})\/(\d{2})\/(\d{4})/;@m2=$b=~m/(\d{2})\/(\d{2})\/(\d{4})/;return -1 if ($m1[2]<$m2[2]);return 1 if ($m1[2]>$m2[2]);return -1 if ($m1[1]<$m2[1]);return 1 if ($m1[1]>$m2[1]);return -1 if ($m1[0]<$m2[0]);return 1 if ($m1[0]>$m2[0]);return 0;}}{@a=split(/\t/);$id=shift(@a);@b=sort compare @a;print $id,$b[0];}' > "$ctl_excl_min" # no header
    "$filter_script" "$ctl_excl_min" "$visit_date" "$instance" > "$ctl_excl_final_icd"
    # join with visit dates, compare two dates, output sample ID if the first date <= second date
    # join -1 1 -2 1 -t$'\t' -o 1.1,1.2,2.2 <(sort -k1,1 "$ctl_excl_min") <(sort -k1,1 "$visit_date") | perl -lne '@a=split(/\t/);$id=shift(@a);@m1=$a[0]=~m/(\d{2})\/(\d{2})\/(\d{4})/;@m2=$a[1]=~m/(\d{2})\/(\d{2})\/(\d{4})/;next if ($m1[2]>$m2[2]);if ($m1[2]<$m2[2]){print $id;next;}next if ($m1[1]>$m2[1]);if ($m1[1]<$m2[1]){print $id;next;}next if ($m1[0]>$m2[0]);if ($m1[0]<=$m2[0]){print $id;}' > "$ctl_excl_final_icd"
else
    : > "$ctl_excl_final_icd"
fi

#------------------------------------------------------------------------------------------------------------------
# 09
# control exclusion OPCS4 codes

ctl_excl_opcs4="$tmpdir"/09.1.ctl_exclusion_opcs4 # with header
ctl_excl_final_opcs4="$tmpdir"/09.2.ctl_exclusion_final_opcs4 # no header

echo $(date "+%d-%b-%Y:%H-%M-%S:") "09 excluding samples from controls based on OPCS4 codes; saving results in $ctl_excl_final_opcs4"

"$date_op_script" --opcs4 "$op_control_exclusion_file" -p "OA" -r "$hesin_release" -o "$ctl_excl_opcs4" >>"$logfile" 2>>"$logfile"

if [[ -e "$ctl_excl_opcs4" ]];then
    tail -n +2  "$ctl_excl_opcs4" | sponge  "$ctl_excl_opcs4"
    "$filter_script" "$ctl_excl_opcs4" "$visit_date" "$instance" > "$ctl_excl_final_opcs4"
    # join with visit dates, compare two dates, output sample ID if the first date <= second date
    # join -1 1 -2 1 -t$'\t' -o 1.1,1.2,2.2 <(tail -n +2 "$ctl_excl_opcs4" | sort -k1,1) <(sort -k1,1 "$visit_date") | perl -lne '@a=split(/\t/);$id=shift(@a);@m1=$a[0]=~m/(\d{2})\/(\d{2})\/(\d{4})/;@m2=$a[1]=~m/(\d{2})\/(\d{2})\/(\d{4})/;next if ($m1[2]>$m2[2]);if ($m1[2]<$m2[2]){print $id;next;}next if ($m1[1]>$m2[1]);if ($m1[1]<$m2[1]){print $id;next;}next if ($m1[0]>$m2[0]);if ($m1[0]<=$m2[0]){print $id;}' > "$ctl_excl_final_opcs4"
else
    : > "$ctl_excl_final_opcs4"
fi

#------------------------------------------------------------------------------------------------------------------
# 10
# final controls

ctl_excl_final="$tmpdir"/10.1.ctl_exclusion_final # no header
ctl_final="$tmpdir"/10.2.control_final # no header

echo $(date "+%d-%b-%Y:%H-%M-%S:") "10 creating final set of controls; saving results in $ctl_final"

cat "$oa_main_out" "$ctl_excl_final_icd" "$ctl_excl_final_opcs4" | sort | uniq > "$ctl_excl_final"
n=$(cat "$ctl_excl_final" | wc -l)
if [[ $n -ne 0 ]];then
    tail -n +2 "$total_samples" | perl -snle 'BEGIN{%h=();if (length($f)!=0){open(fh,"<",$f);while(<fh>){chomp;$h{$_}=1;}close(fh);}}{@a=split(/\t/);if (!defined($h{$a[0]})){print $_;}}' -- -f="$ctl_excl_final" > "$ctl_final"
    # join -t $'\t' -1 1 -2 1 -a 1 -a 2 -e "NA" -o 1.1,2.1 <(tail -n +2 "$total_samples" < sort) <(sort "$ctl_excl_final") | awk 'BEGIN{FS=OFS="\t";}{if ($2=="NA"){print $1;}}' > "$ctl_final"
else
    tail -n +2 "$total_samples" > "$ctl_final"
fi

#------------------------------------------------------------------------------------------------------------------
# 11
# output

echo $(date "+%d-%b-%Y:%H-%M-%S:") "11 creating output; saving results in $outfile"

n=$(cat "$final_cases" | wc -l)
m=$(cat "$ctl_final" | wc -l)
cat <(paste "$final_cases" <(yes "1" | head -n "$n")) <(paste "$ctl_final" <(yes "0" | head -n "$m")) > "$outfile"

#------------------------------------------------------------------------------------------------------------------

if [[ "$keep" == "NO" ]];then
    echo $(date "+%d-%b-%Y:%H-%M-%S:") "deleting temp dir $tmpdir"
    rm -rf "$tmpdir"
fi

date "+%d-%b-%Y:%H-%M-%S" | tee -a "$logfile"

exit 0



