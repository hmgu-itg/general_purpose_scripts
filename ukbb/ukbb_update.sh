#!/usr/bin/env bash

# ALL FILES ARE TAB SEPARATED

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

function usage () {
    echo ""
    echo "Updating script for UKBB data"
    echo ""
    echo "Usage: ukbb_update.sh -f <optional: ID field name; default: \"f.eid\">"
    echo "                      -i input.tab"
    echo "                      -u update.tab"
    echo "                      -o <output dir>"
    echo "                      -b <optional: basename of the output file; default: \"phenotypes\">"
    echo "                      -x <optional: list of individual IDs to exclude>"
    echo "                      -d <optional: debug mode: do not remove temporary files>"
    echo "                      -t <optional: temp dir; default: /tmp>"
    echo "                      -s <optional: save intersecting columns; default: false>"
    echo "                      -p <optional: split bigger file into smaller parts; default: 1>"
    echo "                      -k <optional: for common columns, keep old values for samples not present in the update file; default: false, i.e., report \"NA\">"
    echo ""
    echo "This script updates the input file, which must contain RELEASE/CREATED columns"
    echo "Output release equals input release+1"
    echo "All input/output files are tab-separated"
    echo "Output file contains one header line, with additional CREATED and RELEASE columns"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

declare -A tmp_cats
declare -a tmp_fnames
declare -a tmp_idcol

infile=""
ufile=""
id_field="f.eid"
datestr=$(date +%d-%b-%Y)
out_dir=""
exclude_list=""
debug="NO"
bname="phenotypes"
sorttemp="/tmp"
savex="NO"
keepold="NO"
xname=""
parts=1
while getopts "hi:u:f:o:x:db:t:sp:k" opt; do
    case $opt in
        i)infile=($OPTARG);;
        u)ufile=($OPTARG);;
        f)id_field=($OPTARG);;
        o)out_dir=($OPTARG);;
        x)exclude_list=($OPTARG);;
        d)debug="YES";;
        b)bname=($OPTARG);;
        t)sorttemp=($OPTARG);;
        s)savex="YES";;
        p)parts=($OPTARG);;
        k)keepold="YES";;
        h)usage;;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

exitIfEmpty "$infile" "ERROR: no input file specified"
exitIfEmpty "$ufile" "ERROR: no update file specified"
exitIfEmpty "$out_dir" "ERROR: no output dir specified"
exitIfNotDir "$out_dir" "ERROR: output dir $out_dir is not a directory"

if [[ ! -w "$out_dir" ]];then
    echo "ERROR: output dir $out_dir is not writable" 1>&2
    exit 1
fi

if [[ ! -z "$exclude_list" ]];then
    exitIfNotFile "$exclude_list" "ERROR: exclude list $exclude_list is not a file"
fi

out_dir=${out_dir%/}
sorttemp=${sorttemp%/}

icat=$(getCatCmd "${infile}")
ucat=$(getCatCmd "${ufile}")
release=$(getRelease "${infile}" "${icat}")
if [[ -z "$release" ]];then
    echo "ERROR: could not find RELEASE in $infile" 1>&2
    exit 1
else
    release=$((release+1))
fi

outfile="${out_dir}/${bname}_r${release}.txt.gz"
exitIfExists "$outfile" "ERROR: output file $outfile already exists"
logfile="${out_dir}/${bname}_r${release}.log"

if [[ "$savex" == "YES" ]];then
    xname="${out_dir}/${bname}_r${release}.intersection.txt.gz"
fi

: > "$logfile"

#---------------------------------- CHECKS --------------------------------------

# checking if all rows have the same number of fields
checkFields "$infile" "$icat" "$logfile"
checkDuplicatesHeader "$infile" "$icat" "$logfile"
checkRow "$infile" 1 "$icat" "$logfile"
echo "" | tee -a "$logfile"

checkFields "$ufile" "$ucat" "$logfile"
checkDuplicatesHeader "$ufile" "$ucat" "$logfile"
checkRow "$ufile" 1 "$ucat" "$logfile"
echo "" | tee -a "$logfile"

#----------------------------------------------------

# get ID column number
input_ID_column=$(getColNum "$infile" "$id_field" "$icat")
exitIfEmpty "$input_ID_column" "ERROR: no \"$id_field\" found in $infile"
input_nrows=$("$icat" "$infile"|wc -l)
input_ncols=$("$icat" "$infile"|head -n 1|tr '\t' '\n'|wc -l)
checkColumn "$infile" $input_ID_column $icat "$logfile"
checkDuplicatesColumn "$infile" $input_ID_column $icat "$logfile"
echo ""  |tee -a "$logfile"

update_ID_column=$(getColNum "$ufile" "$id_field" "$ucat")
exitIfEmpty "$update_ID_column" "ERROR: no \"$id_field\" found in $ufile"
update_nrows=$("$ucat" "$ufile"|wc -l)
update_ncols=$("$ucat" "$ufile"|head -n 1|tr '\t' '\n'|wc -l)
checkColumn "$ufile" $update_ID_column $ucat "$logfile"
checkDuplicatesColumn "$ufile" $update_ID_column $ucat "$logfile"
echo ""  |tee -a "$logfile"

#----------------------------------------------------

tmp_fnames+=("$infile")
tmp_fnames+=("$ufile")
tmp_cats["$infile"]="$icat"
tmp_cats["$ufile"]="$ucat"
tmp_idcol+=("$input_ID_column")
tmp_idcol+=("$update_ID_column")

report_missing tmp_fnames tmp_cats tmp_idcol "$logfile"

# REPORT
echo -e "\n---------------------------------------------------------\n" | tee -a "$logfile"

date "+%F %H-%M-%S" | tee -a "$logfile"
echo "Current dir: ${PWD}" | tee -a "$logfile"
echo "Command line: $scriptname ${args[@]}" | tee -a "$logfile"
echo "" | tee -a "$logfile"
echo "INFO: input file: $infile" | tee -a "$logfile"
echo "INFO: rows: $input_nrows" | tee -a "$logfile"
echo "INFO: columns: $input_ncols" | tee -a "$logfile"
echo "" | tee -a "$logfile"
echo "INFO: update file: $ufile" | tee -a "$logfile"
echo "INFO: rows: $update_nrows" | tee -a "$logfile"
echo "INFO: columns: $update_ncols" | tee -a "$logfile"
echo "" | tee -a "$logfile"
echo "OUTPUT RELEASE: $release" | tee -a "$logfile"
echo "KEEP TEMP FILES: $debug" | tee -a "$logfile"
echo "SAVE COMMON COLUMNS: $savex" | tee -a "$logfile"
echo "KEEP OLD VALUES FOR COMMON COLUMNS / MISSING SAMPLES IN THE UPDATE FILE: $keepold" | tee -a "$logfile"
echo "TEMP DIR: $sorttemp" | tee -a "$logfile"
echo "ID FIELD: $id_field" | tee -a "$logfile"
echo "EXCLUDE LIST: $exclude_list" | tee -a "$logfile"
echo "OUTPUT FILE: $outfile" | tee -a "$logfile"
echo "LOG FILE: $logfile" | tee -a "$logfile"

echo -e "\n---------------------------------------------------------\n" | tee -a "$logfile"

#-------------------------------------- CREATING OUTPUT -------------------------------------------------

# FULL OUTER JOIN INPUT FILES ON IDs, THEN EXCLUDE IDS FROM -x LIST 

update_file "$infile" "$ufile" "$sorttemp" "$logfile" tmpf "$xname" "$keepold" "$parts"
if [[ ! -z "$tmpf" ]];then
    # header line
    paste <(head -n 1 "$tmpf") <(echo RELEASE CREATED | tr ' ' '\t') | gzip - -c > "${outfile}"
    # adding output body, excluding IDs from the exclude list
    x=$(cat $tmpf|wc -l)
    x=$((x-1)) # lines in tmpf body
    paste <(tail -n +2 "$tmpf") <(yes $release $datestr | tr ' ' '\t' | head -n $x) | perl -snle 'BEGIN{%h=();if (length($f)!=0){open(fh,"<",$f);while(<fh>){chomp;$h{$_}=1;}close(fh);}}{@a=split(/\t/);if (!defined($h{$a[0]})){print $_;}}' -- -f="$exclude_list" | gzip - -c >> "${outfile}"
    if [[ $debug == "NO" ]];then
	rm $tmpf
    fi
    
    final_rows=$(zcat "${outfile}" | wc -l)
    final_cols=$(zcat "${outfile}" | head -n 1 | tr '\t' '\n' | wc -l)
    echo "INFO: rows in output: $final_rows (including header row)" | tee -a "$logfile"
    echo "INFO: columns in output: $final_cols (including ID column and CREATED, RELEASE columns)" | tee -a "$logfile"
    echo "" | tee -a "$logfile"    
fi

date "+%F %H-%M-%S"  | tee -a "$logfile"
    
exit 0

