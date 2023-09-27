#!/usr/bin/env bash

# ALL FILES ARE TAB SEPARATED

scriptname=$0
args=("$@")

# help formatting
bold=$(tput bold)
underlined=$(tput smul)
normal=$(tput sgr0)

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

function merge_two_files {
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

    local flag
    local gflag=0
    local cat1=$(getCatCmd "$fname1")
    local cat2=$(getCatCmd "$fname2")
    local idCol1=$(getColNum "$fname1" "f.eid" "$cat1")
    local idCol2=$(getColNum "$fname2" "f.eid" "$cat2")
    local ncols1=$("$cat1" "$fname1" | head -n 1 | tr '\t' '\n' | wc -l)
    local ncols2=$("$cat2" "$fname2" | head -n 1 | tr '\t' '\n' | wc -l)

    getColNumbers "$fname1" "$cat1" colnames1
    getColNumbers "$fname2" "$cat2" colnames2

    echo "INFO: merging"  | tee -a "$logfile"
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
    
    if [[ "${#common_cols[@]}" -eq 0 ]];then # no common colnames
	awk 'BEGIN{FS=OFS="\t";}{if ($2=="NA"){$2=$1;}print $0;}' "$tmpfile" | cut -f 2- | sponge "$tmpfile"
	ret="$tmpfile"
    else # there are common colnames
	exclude_cols=(2)
	for c in "${common_cols[@]}";do
	    getColNums "$tmpfile" "$c" "cat" tmp_ar
	    echo -n "INFO: checking common column $c (" $(join_by "," "${tmp_ar[@]}") ") ... " | tee -a "$logfile"
	    flag=$(cut -f 1,2,"${tmp_ar[0]}","${tmp_ar[1]}" "$tmpfile" | awk 'BEGIN{FS="\t";f="OK";}{if ($1!="NA" && $2!="NA" && $3!="NA" && $4!="NA" && $3!=$4){f=$1" "$3" "$4;exit;}}END{print f;}')
	    if [[ "$flag" != "OK" ]];then
		gflag=1
		echo "\nERROR: conflicting values for column $c: $flag\n" | tee -a "$logfile"
	    else
		echo "OK" | tee -a "$logfile"
		awk -v i="${tmp_ar[0]}" -v j="${tmp_ar[1]}" 'BEGIN{FS="\t";}{if ($i=="NA"){$i=$j;}print $0;}' "$tmpfile" | sponge "$tmpfile"
		exclude_cols+=("${tmp_ar[1]}")
	    fi	    
	done
	if [[ "$gflag" -eq 0 ]];then
	    awk 'BEGIN{FS=OFS="\t";}{if ($2=="NA"){$2=$1;}print $0;}' "$tmpfile" | sponge "$tmpfile"
	    s=$(join_by "," "${exclude_cols[@]}")
	    cut --complement -f "$s" "$tmpfile" | sponge "$tmpfile"
	    ret="$tmpfile"
	else
	    rm "$tmpfile"
	    ret=""
	fi
    fi
}


function usage () {
    echo ""
    echo "Merging script for UKBB data"
    echo ""
    echo "Usage: ukbb_merge.sh -f <ID field name; default: \"f.eid\">"
    echo "                     -i input1.tab"
    echo "                     -i input2.tab"
    echo "           ...                    "
    echo "                     -i inputN.tab"
    echo ""
    echo "                     -o <output dir>"
    echo "                     -b <basename of the output file; default: \"phenotypes\">"
    echo "                     -r <release; default: \"1\">"
    echo "                     -x <list of individual IDs to exclude>"
    echo "                     -d <debug mode: do not remove temporary files>"
    echo "                     -a <add CREATED and RELEASE columns to output; default: false>"
    echo "                     -t <temp dir; default: /tmp>"
    echo ""
    echo "This script merges all input files"
    echo "All input/output files are tab-separated"
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
date "+%F %H-%M-%S" | tee -a "$logfile"
echo "Current dir: ${PWD}" | tee -a "$logfile"
echo "Command line: $scriptname ${args[@]}" | tee -a "$logfile"
echo "" | tee -a "$logfile"
echo "Output release: $release" | tee -a "$logfile"
echo "Output file: $outfile" | tee -a "$logfile"
echo "Exclude list: $exclude_list" | tee -a "$logfile"
echo "" | tee -a "$logfile"

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
    while [[ ! -z "$tmpf" && $i -lt $n_input ]];do
	merge_two_files "$tmpf" "${input_fnames[$i]}" "$sorttemp" "$logfile" tmpf
	tmpfiles+=("$tmpf")
	i=$((i+1))
    done

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
	    rm $f
	done	
    fi
fi

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

