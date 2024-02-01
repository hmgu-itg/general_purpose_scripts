#!/usr/bin/env bash

# ALL FILES ARE TAB SEPARATED
# NO COMMON FIELDS ALLOWED

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

function usage () {
    echo ""
    echo "Merging script for UKBB data"
    echo ""
    echo "Usage: ukbb_merge3.sh -f <optional: ID field name; default: \"f.eid\">"
    echo "                      -i input1.tab"
    echo "                      -i input2.tab"
    echo "           ...                    "
    echo "                      -i inputN.tab"
    echo ""
    echo "                      -o <output dir>"
    echo "                      -b <optional: basename of the output file; default: \"phenotypes\">"
    echo "                      -r <optional: release; default: \"1\">"
    echo "                      -x <optional: list of individual IDs to exclude>"
    echo "                      -d <optional: debug mode: do not remove temporary files; default: false>"
    echo "                      -a <optional: add CREATED and RELEASE columns to output; default: false>"
    echo "                      -t <optional: temp dir; default: /tmp>"
    echo ""
    echo "This script merges all input files"
    echo "All input/output files are tab-separated"
    echo "Fields in input files must be disjoint (except for the sample ID field)"
    echo "Output file contains one header line (with additional CREATED and RELEASE columns if \"-a\" is given)"
    echo ""
    exit 0
}

declare -a input_fnames
declare -a input_ID_column
declare -a input_nrows
declare -a input_ncols
declare -A cats
declare -a tmpfiles

if [[ $# -eq 0 ]];then
    usage
fi

id_field="f.eid"
datestr=$(date +%d-%b-%Y)
out_dir=""
release=""
exclude_list=""
debug="NO"
add_release="NO"
bname="phenotypes"
sorttemp="/tmp"
while getopts "hi:f:o:r:x:dab:t:" opt; do
    case $opt in
        i)input_fnames+=($OPTARG);;
        f)id_field=($OPTARG);;
        o)out_dir=($OPTARG);;
        r)release=($OPTARG);;
        x)exclude_list=($OPTARG);;
        d)debug="YES";;
        a)add_release="YES";;
        b)bname=($OPTARG);;
        t)sorttemp=($OPTARG);;
        h)usage;;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

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

n_input="${#input_fnames[@]}"

if [[ $n_input -eq 0 ]];then
    echo "ERROR: no input files specified" 1>&2
    exit 1
fi

# get zcat/cat command for each input file
for i in $(seq 0 $((n_input-1)));do
    cats["${input_fnames[$i]}"]=$(getCatCmd "${input_fnames[$i]}")
done

# if no release specified, set it to "1"
if [[ -z "$release" ]];then
    release="1"
fi

outfile="${out_dir}/${bname}_r${release}.txt.gz"
exitIfExists "$outfile" "ERROR: output file $outfile already exists"
logfile="${out_dir}/${bname}_r${release}.log"

: > "$logfile"

#----------------------------------------------------

# checking if all rows have the same number of fields
for i in $(seq 0 $((n_input-1)));do
    checkFields "${input_fnames[$i]}" "${cats[${input_fnames[$i]}]}" "$logfile"
    checkDuplicatesHeader "${input_fnames[$i]}" "${cats[${input_fnames[$i]}]}" "$logfile"
    checkRow "${input_fnames[$i]}" 1 "${cats[${input_fnames[$i]}]}" "$logfile"
    echo "" | tee -a "$logfile"
done

#----------------------------------------------------

# get ID column number for every input file
for i in $(seq 0 $((n_input-1)));do
    x=$(getColNum "${input_fnames[$i]}" "$id_field" ${cats[${input_fnames[$i]}]})
    exitIfEmpty "$x" "ERROR: no \"$id_field\" found in ${input_fnames[$i]}"
    input_ID_column+=($x)
    input_nrows+=($(${cats[${input_fnames[$i]}]} "${input_fnames[$i]}"|wc -l))
    input_ncols+=($(${cats[${input_fnames[$i]}]} "${input_fnames[$i]}"|head -n 1|tr '\t' '\n'|wc -l))
    checkColumn "${input_fnames[$i]}" $input_ID_column ${cats[${input_fnames[$i]}]} "$logfile"
    checkDuplicatesColumn "${input_fnames[$i]}" $input_ID_column ${cats[${input_fnames[$i]}]} "$logfile"
    echo ""  |tee -a "$logfile"
done

#----------------------------------------------------

# REPORT
for i in $(seq 0 $((n_input-1)));do
    echo "INFO: input file: ${input_fnames[$i]}" | tee -a "$logfile"
    echo "INFO: rows: ${input_nrows[$i]}" | tee -a "$logfile"
    echo "INFO: columns: ${input_ncols[$i]}" | tee -a "$logfile"
    echo "" | tee -a "$logfile"
done

report_missing input_fnames cats input_ID_column "$logfile"

#----------------------------------------------------

#
# check if column names (except for the ID field) in input files are unique
#

echo "Checking if column names in input files are disjoint ... " | tee -a "$logfile"
flag=0
if [[ $n_input -gt 1 ]];then
    for i in $(seq 0 $((n_input-1)));do
	for j in $(seq $((i+1)) $((n_input-1)));do
	    x=$(cat <("${cats[${input_fnames[$i]}]}" "${input_fnames[$i]}" | head -n 1 | cut --complement -f ${input_ID_column[$i]} | tr '\t' '\n') <("${cats[${input_fnames[$j]}]}" "${input_fnames[$j]}" | head -n 1 | cut --complement -f ${input_ID_column[$j]} | tr '\t' '\n') | sort -T "${sorttemp}" | uniq -d | wc -l)
	    if [[ $x -ne 0 ]];then
		flag=1
		echo "ERROR: input files ${input_fnames[$i]} and ${input_fnames[$j]} have columns in common" | tee -a "$logfile"
	    else
		echo "INFO: files $i $j: OK" | tee -a "$logfile"
	    fi	    
	done
    done
fi

if [[ "$flag" -ne 0 ]];then
    exit 1
fi

echo -e "\n---------------------------------------------------------\n" | tee -a "$logfile"

date "+%F %H-%M-%S" | tee -a "$logfile"
echo "Current dir: ${PWD}" | tee -a "$logfile"
echo "Command line: $scriptname ${args[@]}" | tee -a "$logfile"
echo "" | tee -a "$logfile"
echo "INPUT FILES: $n_input" | tee -a "$logfile"
for i in $(seq 0 $((n_input-1)));do
    echo "${input_fnames[$i]}" | tee -a "$logfile"
done
echo "" | tee -a "$logfile"
echo "OUTPUT RELEASE: $release" | tee -a "$logfile"
echo "ADD RELEASE: $add_release" | tee -a "$logfile"
echo "KEEP TEMP FILES: $debug" | tee -a "$logfile"
echo "TEMP DIR: $sorttemp" | tee -a "$logfile"
echo "ID FIELD: $id_field" | tee -a "$logfile"
echo "EXCLUDE LIST: $exclude_list" | tee -a "$logfile"
echo "OUTPUT FILE: $outfile" | tee -a "$logfile"
echo "LOG FILE: $logfile" | tee -a "$logfile"

echo -e "\n---------------------------------------------------------\n" | tee -a "$logfile"

#-------------------------------------- CREATING OUTPUT -------------------------------------------------

# OUTER JOIN INPUT FILES ON IDs, THEN EXCLUDE IDS FROM -x LIST 

if [[ $n_input -eq 1 ]];then # only one input file
    if [[ "$add_release" == "YES" ]];then
	paste <(${cats["${input_fnames[0]}"]} "${input_fnames[0]}" | head -n 1) <(echo RELEASE CREATED|tr ' ' '\t') | gzip - -c > "${outfile}"
	${cats["${input_fnames[0]}"]} "${input_fnames[0]}" | tail -n +2 | perl -snle 'BEGIN{$,="\t";%h=();if (length($f)!=0){open(fh,"<",$f);while(<fh>){chomp;$h{$_}=1;}close(fh);}}{@a=split(/\t/);if (!defined($h{$a[$c-1]})){print $_,$r,$d;}}' -- -f="$exclude_list" -c="${input_ID_column[0]}" -r=$release -d=$datestr| gzip - -c >> "${outfile}"
    else # no release info in output
	${cats["${input_fnames[0]}"]} "${input_fnames[0]}" | head -n 1 |gzip - -c > "${outfile}"
	${cats["${input_fnames[0]}"]} "${input_fnames[0]}" | tail -n +2 | perl -snle 'BEGIN{$,="\t";%h=();if (length($f)!=0){open(fh,"<",$f);while(<fh>){chomp;$h{$_}=1;}close(fh);}}{@a=split(/\t/);if (!defined($h{$a[$c-1]})){print $_;}}' -- -f="$exclude_list" -c="${input_ID_column[0]}" | gzip - -c >> "${outfile}"
    fi
else # several input files
    tmpfiles=()
    # OUTER JOIN on IDs
    merge_two_files "${input_fnames[0]}" "${input_fnames[1]}" "$sorttemp" "$logfile" tmpf
    tmpfiles+=("$tmpf")
    i=2
    while [[ ( ! -z "$tmpf" ) && ( $i -lt $n_input ) ]];do
	merge_two_files "$tmpf" "${input_fnames[$i]}" "$sorttemp" "$logfile" tmpf
	tmpfiles+=("$tmpf")
	i=$((i+1))
	# echo "DEBUG: i=$i"
	# echo "DEBUG: tmpf=$tmpf"
    done
    echo "INFO: done merging" | tee -a "$logfile"

    if [[ ! -z "$tmpf" ]];then
	# header line
	if [[ "$add_release" == "YES" ]];then
	    paste <(head -n 1 "$tmpf") <(echo RELEASE CREATED | tr ' ' '\t') | gzip - -c > "${outfile}"
	else
	    head -n 1 "$tmpf" | gzip - -c > "${outfile}"
	fi
	# adding output body, excluding IDs from the exclude list
	x=$(cat $tmpf|wc -l)
	x=$((x-1)) # lines in tmpf body
	if [[ "$add_release" == "YES" ]];then
	    paste <(tail -n +2 "$tmpf") <(yes $release $datestr | tr ' ' '\t' | head -n $x) | perl -snle 'BEGIN{%h=();if (length($f)!=0){open(fh,"<",$f);while(<fh>){chomp;$h{$_}=1;}close(fh);}}{@a=split(/\t/);if (!defined($h{$a[0]})){print $_;}}' -- -f="$exclude_list" | gzip - -c >> "${outfile}"
	else
	    tail -n +2 "$tmpf" | perl -snle 'BEGIN{%h=();if (length($f)!=0){open(fh,"<",$f);while(<fh>){chomp;$h{$_}=1;}close(fh);}}{@a=split(/\t/);if (!defined($h{$a[0]})){print $_;}}' -- -f="$exclude_list" | gzip - -c >> "${outfile}"
	fi
    fi

    if [[ $debug == "NO" ]];then
	for f in "${tmpfiles[@]}";do
	    echo "INFO: deleting temporary file: $f" | tee -a "$logfile"
	    rm $f
	done	
    fi
fi

echo "" | tee -a "$logfile"    

final_rows=$(zcat "${outfile}" | wc -l)
final_cols=$(zcat "${outfile}" | head -n 1 | tr '\t' '\n' | wc -l)
echo "INFO: rows in output: $final_rows (including header row)" | tee -a "$logfile"
    if [[ "$add_release" == "YES" ]];then
	echo "INFO: columns in output: $final_cols (including ID column and CREATED, RELEASE columns)" | tee -a "$logfile"
    else
	echo "INFO: columns in output: $final_cols" | tee -a "$logfile"
    fi
echo "" | tee -a "$logfile"    

date "+%F %H-%M-%S"  | tee -a "$logfile"
    
exit 0

