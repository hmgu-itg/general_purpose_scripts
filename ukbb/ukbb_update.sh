#!/usr/bin/env bash

# ALL FILES ARE TAB SEPARATED

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

function update_file {
    local fname1=$1
    local fname2=$2
    local tmpdir=$3
    local logfile=$4
    local -n ret=$5

    declare -A colnames1
    declare -A colnames2
    declare -a common_cols
    declare -a tmp_ar
    declare -a exclude_cols

    local cat1=$(getCatCmd "$fname1")
    local cat2=$(getCatCmd "$fname2")
    local idCol1=$(getColNum "$fname1" "f.eid" "$cat1")
    local idCol2=$(getColNum "$fname2" "f.eid" "$cat2")
    local ncols1=$("$cat1" "$fname1" | head -n 1 | tr '\t' '\n' | wc -l)
    local ncols2=$("$cat2" "$fname2" | head -n 1 | tr '\t' '\n' | wc -l)

    getColNumbers "$fname1" "$cat1" colnames1
    getColNumbers "$fname2" "$cat2" colnames2

    echo "INFO: updating"  | tee -a "$logfile"
    echo "INFO: file1: $fname1"  | tee -a "$logfile"
    echo "INFO: file2: $fname2"  | tee -a "$logfile"
    echo ""  | tee -a "$logfile"
    
    for c in "${!colnames1[@]}";do
	if [[ -v "colnames2[$c]" && "$c" != "f.eid" ]];then
	    common_cols+=($c)
	fi
    done

    echo "INFO: common fields: ${#common_cols[@]}"  | tee -a "$logfile"
    echo ""  | tee -a "$logfile"
    
    tmpfile=$(mktemp -p "$tmpdir" merge_XXXXXXXX)
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not create tmpfile" | tee -a "$logfile"
	ret=""
	return
    fi

    echo "INFO: output: $tmpfile"  | tee -a "$logfile"
    echo ""  | tee -a "$logfile"
    
    fmt="1.${idCol1},2.${idCol2}"
    for i in $(seq 2 ${ncols1});do
	fmt=$fmt",1.$i"
    done
    for i in $(seq 2 ${ncols2});do
	fmt=$fmt",2.$i"
    done

    echo "INFO: joining input files"  | tee -a "$logfile"
    echo ""  | tee -a "$logfile"
    
    join_cmd="join --header -t$'\t' -1 ${idCol1} -2 ${idCol2} -a 1 -a 2 -e NA -o $fmt <(cat <(${cat1} ${fname1} | head -n 1) <(${cat1} ${fname1} | tail -n +2 | sort -T ${tmpdir} -t$'\t' -k${idCol1},${idCol1})) <(cat <(${cat2} ${fname2} | head -n 1) <(${cat2} ${fname2} | tail -n +2 | sort -T ${tmpdir} -t$'\t' -k${idCol2},${idCol2}))"
    eval "$join_cmd > $tmpfile"

    x=$(tail -n +2  "$tmpfile" | awk 'BEGIN{FS="\t";c=0;}{if ($1!="NA" && $2=="NA"){c=c+1;}}END{print c;}')
    echo "INFO: samples in file1 only: $x"  | tee -a "$logfile"
    x=$(tail -n +2  "$tmpfile" | awk 'BEGIN{FS="\t";c=0;}{if ($1=="NA" && $2!="NA"){c=c+1;}}END{print c;}')
    echo "INFO: samples in file2 only: $x"  | tee -a "$logfile"
    x=$(tail -n +2  "$tmpfile" | awk 'BEGIN{FS="\t";c=0;}{if ($1!="NA" && $2!="NA"){c=c+1;}}END{print c;}')
    echo "INFO: samples in both file1 and file2: $x"  | tee -a "$logfile"
    x=$(head -n 1 "$tmpfile" | tr '\t' '\n' | wc -l )
    echo "INFO: total columns in joined file: $x"  | tee -a "$logfile"
    echo ""  | tee -a "$logfile"

    echo "DEBUG: excluding second column"  | tee -a "$logfile"
    awk 'BEGIN{FS=OFS="\t";}{if ($2=="NA"){$2=$1;}print $0;}' "$tmpfile" | cut -f 2- | TMPDIR="${tmpdir}" sponge "$tmpfile"
    echo "DEBUG: done"  | tee -a "$logfile"
    # tmpfile contains one ID column (1st), all columns from file1 (including CREATED, RELEASE), all columns from file2
    
    if [[ "${#common_cols[@]}" -eq 0 ]];then # no common colnames, remove RELEASE, CREATED columns
	cut --complement -f $(getColNum "$tmpfile" "RELEASE"),$(getColNum "$tmpfile" "CREATED") "$tmpfile" | TMPDIR="${tmpdir}" sponge "$tmpfile"
	ret="$tmpfile"
    else # there are common colnames
	# report some stats comparing shared columns in both files
	"${scriptdir}/get_stats.pl" "$tmpfile" >> "$logfile"
	# remove RELEASE, CREATED columns
	exclude_cols=()
	exclude_cols+=($(getColNum "$tmpfile" "RELEASE"))
	exclude_cols+=($(getColNum "$tmpfile" "CREATED"))
	# remove common columns from the first file
	for c in "${common_cols[@]}";do
	    getColNums "$tmpfile" "$c" "cat" tmp_ar
	    exclude_cols+=("${tmp_ar[0]}")
	done
	cut --complement -f $(join_by "," "${exclude_cols[@]}") "$tmpfile" | TMPDIR="${tmpdir}" sponge "$tmpfile"
	ret="$tmpfile"
    fi
}

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

infile=""
ufile=""
id_field="f.eid"
datestr=$(date +%d-%b-%Y)
out_dir=""
exclude_list=""
debug="NO"
bname="phenotypes"
sorttemp="/tmp"
while getopts "hi:u:f:o:x:db:t:" opt; do
    case $opt in
        i)infile=($OPTARG);;
        u)ufile=($OPTARG);;
        f)id_field=($OPTARG);;
        o)out_dir=($OPTARG);;
        x)exclude_list=($OPTARG);;
        d)debug="YES";;
        b)bname=($OPTARG);;
        t)sorttemp=($OPTARG);;
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

: > "$logfile"
date "+%F %H-%M-%S" | tee -a "$logfile"
echo "Current dir: ${PWD}" | tee -a "$logfile"
echo "Command line: $scriptname ${args[@]}" | tee -a "$logfile"
echo "" | tee -a "$logfile"
echo "Output release: $release" | tee -a "$logfile"
echo "Output file: $outfile" | tee -a "$logfile"
echo "Exclude list: $exclude_list" | tee -a "$logfile"
echo "" | tee -a "$logfile"

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

# REPORT
echo "INFO: input file: $infile" | tee -a "$logfile"
echo "INFO: rows: $input_nrows" | tee -a "$logfile"
echo "INFO: columns: $input_ncols" | tee -a "$logfile"
echo "" | tee -a "$logfile"

echo "INFO: update file: $ufile" | tee -a "$logfile"
echo "INFO: rows: $update_nrows" | tee -a "$logfile"
echo "INFO: columns: $update_ncols" | tee -a "$logfile"
echo "" | tee -a "$logfile"

#-------------------------------------- CREATING OUTPUT -------------------------------------------------

# OUTER JOIN INPUT FILES ON IDs, THEN EXCLUDE IDS FROM -x LIST 

update_file "$infile" "$ufile" "$sorttemp" "$logfile" tmpf
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

