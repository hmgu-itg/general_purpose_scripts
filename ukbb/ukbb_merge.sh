#!/usr/bin/env bash

# ALL INPUT FILES ARE TAB SEPARATED

scriptname=$0
args=("$@")

# help folrmatting
bold=$(tput bold)
underlined=$(tput smul)
normal=$(tput sgr0)

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

function usage () {
    echo ""
    echo "Merging/updating script for UKBB data"
    echo ""
    echo "Usage: ukbb_merge.sh -f <ID field name; default: \"f.eid\">"
    echo "                     -i input1.tab"
    echo "                     -i input2.tab"
    echo "           ...                    "
    echo "                     -i inputN.tab"
    echo ""
    echo "                     -u update1.tab"
    echo "                     -u update2.tab"
    echo "           ...                    "
    echo "                     -u updateM.tab"
    echo ""
    echo "                     -o <output dir>"
    echo "                     -r <release; default: \"1\" if merging, incrementing RELEASE in input if updating>"
    echo "                     -k <when updating, also include samples from update file(s) that are ${bold}not${normal} present in input; default: false>"
    echo "                     -x <list of IDs to exclude>"
    echo ""
    echo "All input/update/output files are tab-separated"
    echo ""
    echo "${underlined}Merge mode${normal}: if no update (-u) files are specified, the script merges all input files."
    echo ""
    echo "${underlined}Update mode${normal}: if at least one update (-u) file is given, the script only uses the first input (-i) file" 
    echo "and ignores the remaining input files. If a column in the input file"
    echo "is present in an update file, its content in the input file will be updated."
    echo ""
    exit 0
}

declare -a input_fnames
declare -a update_fnames
declare -a input_ID_column
declare -a update_ID_column
declare -a input_nrows
declare -a update_nrows
declare -a input_ncols
declare -a update_ncols
declare -a class_array
declare -A column_class
declare -A cats

if [[ $# -eq 0 ]];then
    usage
fi

id_field="f.eid"
datestr=$(date +%F)
out_dir=""
release=""
exclude_list=""
mode="inner"
while getopts "hi:u:f:o:r:x:k" opt; do
    case $opt in
        i)input_fnames+=($OPTARG);;
        u)update_fnames+=($OPTARG);;
        f)id_field=($OPTARG);;
        o)out_dir=($OPTARG);;
        r)release=($OPTARG);;
        r)exclude_list=($OPTARG);;
        k)mode="right";;
        h)usage;;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

ret=$(checkInstalledCommands)
if [[ $ret -ne 0 ]];then
    exit 1
fi

if [[ -z "$out_dir" ]];then
    echo "ERROR: no output dir specified" 1>&2
    exit 1
fi

if [[ ! -d "$out_dir" ]];then
    echo "ERROR: output dir $out_dir is not a directory" 1>&2
    exit 1
fi

if [[ ! -w "$out_dir" ]];then
    echo "ERROR: output dir $out_dir is not writable" 1>&2
    exit 1
fi

out_dir=${out_dir%/}

n_input="${#input_fnames[@]}"
n_update="${#update_fnames[@]}"

if [[ $n_input -eq 0 ]];then
    echo "ERROR: no input files specified" 1>&2
    exit 1
fi

if [[ $n_input -lt 2 && $n_update -eq 0 ]];then
    echo "ERROR: no update files specified, only one input file specified (need at least two input files to merge)" 1>&2
    exit 1
fi

if [[ $n_input -eq 0 && $n_update -gt 0 ]];then
    echo "ERROR: update file(s) specified but no input files specified" 1>&2
    exit 1
fi

if [[ $n_update -gt 0 && $n_input -gt 1 ]];then
    echo "WARN: $n_update update files and $n_input input files specified; only the first input file (${input_fnames[0]}) will be processed" 1>&2
fi

# get zcat/cat command
for i in $(seq 0 $((n_input-1)));do
    cats["${input_fnames[$i]}"]=$(getCatCmd "${input_fnames[$i]}")
done
for i in $(seq 0 $((n_update-1)));do
    cats["${update_fnames[$i]}"]=$(getCatCmd "${update_fnames[$i]}")
done

# get release, either from command line or from input file, in case of update
if [[ $n_update -eq 0 ]];then
    # if no release specified, set it to "1"
    if [[ -z "$release" ]];then
	release="1"
    fi
else
    # read previous release from input[0]
    release_col=$(getColNum "${input_fnames[0]}" "RELEASE" ${cats[${input_fnames[0]}]})
    if [[ -z "$release_col" ]];then
	echo "ERROR: no RELEASE column found in the input file ${input_fnames[0]}" 1>&2
	exit 1
    fi
    created_col=$(getColNum "${input_fnames[0]}" "CREATED" ${cats[${input_fnames[0]}]})
    if [[ -z "$created_col" ]];then
	echo "ERROR: no CREATED column found in the input file ${input_fnames[0]}" 1>&2
	exit 1
    fi
    if [[ -z "$release" ]];then
	release=$(${cats[${input_fnames[0]}]} "${input_fnames[0]}"|cut -f ${release_col}|head -n 3|tail -n 1)
	release=$((release+1))
    fi
fi

outfile="${out_dir}/phenotypes_r${release}.txt.gz"
if [[ -f "$outfile" ]];then
    echo "ERROR: output file $outfile already exists" 1>&2
    exit 1
fi
logfile="${out_dir}/phenotypes_r${release}.log"

: > "$logfile"
date "+%F %H-%M-%S"|tee -a "$logfile"
echo "Current dir: ${PWD}"|tee -a "$logfile"
echo "Command line: $scriptname ${args[@]}"|tee -a "$logfile"
echo ""|tee -a "$logfile"
echo "Output release: $release"|tee -a "$logfile"
echo "Output file: $outfile"|tee -a "$logfile"
echo "INFO: join mode: $mode"|tee -a "$logfile"
echo ""|tee -a "$logfile"
#----------------------------------------------------

# checking if all rows have the same number of fields
for i in $(seq 0 $((n_input-1)));do
    checkFields "${input_fnames[$i]}" "${cats[${input_fnames[$i]}]}" "$logfile"
    checkDuplicatesHeader "${input_fnames[$i]}" "${cats[${input_fnames[$i]}]}" "$logfile"
    checkRow "${input_fnames[$i]}" 1 "${cats[${input_fnames[$i]}]}" "$logfile"
    echo ""|tee -a "$logfile"
done
for i in $(seq 0 $((n_update-1)));do
    checkFields "${update_fnames[$i]}" "${cats[${update_fnames[$i]}]}" "$logfile"
    checkDuplicatesHeader "${update_fnames[$i]}" "${cats[${update_fnames[$i]}]}" "$logfile"
    checkRow "${update_fnames[$i]}" 1 "${cats[${update_fnames[$i]}]}" "$logfile"
    echo ""|tee -a "$logfile"
done

# check classes row if we're updating
if [[ $n_update -gt 0 ]];then
    checkRow "${input_fnames[0]}" 2 "${cats[${input_fnames[0]}]}" "$logfile"
fi
#----------------------------------------------------

# get ID column 
for i in $(seq 0 $((n_input-1)));do
    x=$(getColNum "${input_fnames[$i]}" "$id_field" ${cats[${input_fnames[$i]}]})
    if [[ -z $x ]];then
	echo "ERROR: no \"$id_field\" found in ${input_fnames[$i]}"|tee -a "$logfile"
	exit 1
    fi    
    input_ID_column+=($x)
    input_nrows+=($(${cats[${input_fnames[$i]}]} "${input_fnames[$i]}"|wc -l))
    input_ncols+=($(${cats[${input_fnames[$i]}]} "${input_fnames[$i]}"|head -n 1|tr '\t' '\n'|wc -l))
    checkColumn "${input_fnames[$i]}" $input_ID_column ${cats[${input_fnames[$i]}]} "$logfile"
    checkDuplicatesColumn "${input_fnames[$i]}" $input_ID_column ${cats[${input_fnames[$i]}]} "$logfile"
    echo ""  |tee -a "$logfile"
done
for i in $(seq 0 $((n_update-1)));do
    x=$(getColNum "${update_fnames[$i]}" "$id_field" ${cats[${update_fnames[$i]}]})
    if [[ -z $x ]];then
	echo "ERROR: no \"$id_field\" found in ${update_fnames[$i]}"|tee -a "$logfile"
	exit 1
    fi    
    update_ID_column+=($x)
    update_nrows+=($(${cats[${update_fnames[$i]}]} "${update_fnames[$i]}"|wc -l))
    update_ncols+=($(${cats[${update_fnames[$i]}]} "${update_fnames[$i]}"|head -n 1|tr '\t' '\n'|wc -l))
    checkColumn "${update_fnames[$i]}" $update_ID_column ${cats[${update_fnames[$i]}]} "$logfile"
    checkDuplicatesColumn "${update_fnames[$i]}" $update_ID_column ${cats[${update_fnames[$i]}]} "$logfile"
    echo ""  |tee -a "$logfile"
done
#----------------------------------------------------

# output info
for i in $(seq 0 $((n_input-1)));do
    echo "INFO: input file: ${input_fnames[$i]}"|tee -a "$logfile"
    echo "INFO: ID column index: ${input_ID_column[$i]}"|tee -a "$logfile"
    echo "INFO: rows: ${input_nrows[$i]}"|tee -a "$logfile"
    echo "INFO: columns: ${input_ncols[$i]}"|tee -a "$logfile"
    echo ""|tee -a "$logfile"
done

for i in $(seq 0 $((n_update-1)));do
    echo "INFO: update file: ${update_fnames[$i]}"|tee -a "$logfile"
    echo "INFO: ID column index: ${update_ID_column[$i]}"|tee -a "$logfile"
    echo "INFO: rows: ${update_nrows[$i]}"|tee -a "$logfile"
    echo "INFO: columns: ${update_ncols[$i]}"|tee -a "$logfile"
    echo ""|tee -a "$logfile"
done
#----------------------------------------------------

#
# check if all input files have same IDs
#
# echo -n "Checking if input files have same IDs ... "|tee -a "$logfile"
# for i in $(seq 1 $((n_input-1)));do
#     x=$(cat <(${cats[${input_fnames[0]}]} "${input_fnames[0]}"|cut -f ${input_ID_column[0]}) <(${cats[${input_fnames[$i]}]} "${input_fnames[$i]}"|cut -f ${input_ID_column[$i]})|sort|uniq -u|wc -l)
#     if [[ $x -ne 0 ]];then
# 	echo "ERROR: files ${input_fnames[0]} and ${input_fnames[$i]} have different sets of IDs"|tee -a "$logfile"
# 	exit 1
#     fi
# done
# echo "OK"|tee -a "$logfile"
#----------------------------------------------------

#
# check if all update files have same IDs
#
# echo -n "Checking if update files have same IDs ... "|tee -a "$logfile"
# for i in $(seq 1 $((n_update-1)));do
#     x=$(cat <("${cats[${update_fnames[0]}]}" "${update_fnames[0]}"|cut -f ${update_ID_column[0]}) <("${cats[${update_fnames[$i]}]}" "${update_fnames[$i]}"|cut -f ${update_ID_column[$i]})|sort|uniq -u|wc -l)
#     if [[ $x -ne 0 ]];then
# 	echo "ERROR: files ${input_fnames[0]} and ${input_fnames[$i]} have different sets of IDs"|tee -a "$logfile"
# 	exit 1
#     fi
# done
# echo "OK"|tee -a "$logfile"
#----------------------------------------------------

#
# check if column names (excepth for ID field) in input files are disjoint
#
echo -n "Checking if column names in input files are disjoint ... "|tee -a "$logfile"
if [[ $n_input -gt 1 ]];then
    for i in $(seq 0 $((n_input-1)));do
	for j in $(seq $((i+1)) $((n_input-1)));do
	    x=$(cat <("${cats[${input_fnames[$i]}]}" "${input_fnames[$i]}"|head -n 1|cut --complement -f ${input_ID_column[$i]}) <("${cats[${input_fnames[$j]}]}" "${input_fnames[$j]}"|head -n 1|cut --complement -f ${input_ID_column[$j]})|sort|uniq -d|wc -l)
	    if [[ $x -ne 0 ]];then
		echo "ERROR: input files ${input_fnames[$i]} and ${input_fnames[$j]} have columns in common"|tee -a "$logfile"
		exit 1
	    fi	    
	done
    done
fi
echo "OK"|tee -a "$logfile"
#----------------------------------------------------

#
# check if column names (excepth for ID field) in update files are disjoint
#
echo -n "Checking if column names in update files are disjoint ... "|tee -a "$logfile"
if [[ $n_update -gt 1 ]];then
    for i in $(seq 0 $((n_update-1)));do
	for j in $(seq $((i+1)) $((n_update-1)));do
	    x=$(cat <("${cats[${update_fnames[$i]}]}" "${update_fnames[$i]}"|head -n 1|cut --complement -f ${update_ID_column[$i]}) <("${cats[${update_fnames[$j]}]}" "${update_fnames[$j]}"|head -n 1|cut --complement -f ${update_ID_column[$j]})|sort|uniq -d|wc -l)
	    if [[ $x -ne 0 ]];then
		echo "ERROR: update files ${update_fnames[$i]} and ${update_fnames[$j]} have columns in common"|tee -a "$logfile"
		exit 1
	    fi	    
	done
    done
fi
echo "OK"|tee -a "$logfile"
echo ""|tee -a "$logfile"

#-------------------------------------- OUTPUT -------------------------------------------------

if [[ $n_update -eq 0 ]];then
    # just merging input files
    
    # column classes of input columns are based on source input files
    for i in $(seq 0 $((n_input-1)));do
    	for c in $("${cats[${input_fnames[$i]}]}" "${input_fnames[$i]}"|head -n 1|cut --complement -f ${input_ID_column[$i]});do
    	    column_class[$c]=$i
    	done
    done
    column_class["$id_field"]="NA"

    cmd="cat <(${cats[${input_fnames[0]}]} ${input_fnames[0]}| tail -n +2|cut -f ${input_ID_column[0]})"
    for i in $(seq 1 $((n_input-1)));do
	cmd=$cmd" <(${cats[${input_fnames[$i]}]} ${input_fnames[$i]}| tail -n +2|cut -f ${input_ID_column[$i]})"
    done
    cmn=$(eval $cmd | sort|uniq -c|sed 's/^  *//'|gawk -v n=$n_input '$1==n{print $0;}'| wc -l)
    echo "INFO: output $cmn common IDs between $n_input input files"|tee -a "$logfile"
    echo ""|tee -a "$logfile"
    echo -n "Merging input files ... "|tee -a "$logfile"
    
    # # HEADER
    # command1="paste <(${cats[${input_fnames[0]}]} ${input_fnames[0]}|head -n 1)"
    # for i in $(seq 1 $((n_input-1)));do
    # 	command1=${command1}" <(${cats[${input_fnames[$i]}]} ${input_fnames[$i]}|head -n 1|cut --complement -f ${input_ID_column[$i]})"
    # done
    # command1=${command1}" <(echo RELEASE CREATED| tr ' ' '\t')"

    # for c in $(${cats[${input_fnames[0]}]} ${input_fnames[0]}|head -n 1);do
    # 	class_array+=(${column_class[$c]})
    # done
    # for i in $(seq 1 $((n_input-1)));do
    # 	for c in $(${cats[${input_fnames[$i]}]} ${input_fnames[$i]}|head -n 1|cut --complement -f ${input_ID_column[$i]});do
    # 	    class_array+=(${column_class[$c]})
    # 	done
    # done
    # class_array+=("NA")
    # class_array+=("NA")
    # header2=$(join_by , ${class_array[*]})
    
    # command2="paste <(${cats[${input_fnames[0]}]} ${input_fnames[0]}|tail -n +2|sort -k${input_ID_column[0]},${input_ID_column[0]})"
    # for i in $(seq 1 $((n_input-1)));do
    # 	command2=${command2}" <(${cats[${input_fnames[$i]}]} ${input_fnames[$i]}|tail -n +2|sort -k${input_ID_column[$i]},${input_ID_column[$i]}|cut --complement -f ${input_ID_column[$i]})"
    # done
    # x=${input_nrows[0]}
    # x=$((x-1))
    # command2=${command2}" <(yes $release $datestr| tr ' ' '\t'| head -n $x)"

    # eval "cat <($command1) <(echo $header2|tr ',' '\t') <($command2) | gzip - -c > ${outfile}"

    # only keeping commong IDs between input files
    tmpfile1=$(mktemp -p "$out_dir" temp_1_join_XXXXXXXX)
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not create tmpfile1 "|tee -a "$logfile"
	exit 1
    fi
    
    join_cmd="join --header -t$'\t' -1 ${input_ID_column[0]} -2 ${input_ID_column[1]} <(cat <(${cats[${input_fnames[0]}]} ${input_fnames[0]}|head -n 1) <(${cats[${input_fnames[0]}]} ${input_fnames[0]}|tail -n +2|sort -k${input_ID_column[0]},${input_ID_column[0]})) <(cat <(${cats[${input_fnames[1]}]} ${input_fnames[1]}|head -n 1) <(${cats[${input_fnames[1]}]} ${input_fnames[1]}|tail -n +2|sort -k${input_ID_column[1]},${input_ID_column[1]}))"
    for i in $(seq 2 $((n_input-1)));do
	join_cmd=$join_cmd"|join --header -t$'\t' -1 1 -2 ${input_ID_column[$i]} - <(cat <(${cats[${input_fnames[$i]}]} ${input_fnames[$i]}|head -n 1) <(${cats[${input_fnames[$i]}]} ${input_fnames[$i]}|tail -n +2|sort -k${input_ID_column[$i]},${input_ID_column[$i]}))"
    done
    eval "$join_cmd > $tmpfile1"

    # adding class line
    paste <(head -n 1 "$tmpfile1") <(echo RELEASE CREATED|tr ' ' '\t')|gzip - -c > "${outfile}"
    for c in $(head -n 1 "$tmpfile1");do
	class_array+=(${column_class[$c]})
    done
    class_array+=("NA" "NA")
    cline=$(join_by , "${class_array[@]}")
    eval "cat <(echo $cline|tr ',' '\t')|gzip - -c >> ${outfile}"    
    x=$(cat $tmpfile1|wc -l)
    x=$((x-1))
    paste <(tail -n +2 "$tmpfile1") <(yes $release $datestr|tr ' ' '\t'|head -n $x)|gzip - -c >> "${outfile}"    
    rm -f "$tmpfile1"
    
    echo "Done"|tee -a "$logfile"
    date "+%F %H-%M-%S" |tee -a "$logfile"
else
    # updating the first input file using update files
    
    echo "Updating input ... "|tee -a "$logfile"

    #----------------------------------------------------------------
    # merging all update files and saving the result in a temporary file, with classes
    # update files have same IDs, so using paste instead of join
    tmpfile0=$(mktemp -p "$out_dir" temp_0_join_XXXXXXXX)
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not create tmpfile0 "|tee -a "$logfile"
	exit 1
    fi
    tmpfile1=$(mktemp -p "$out_dir" temp_1_update_XXXXXXXX)
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not create tmpfile1 "|tee -a "$logfile"
	exit 1
    fi

    cmd="cat <(${cats[${update_fnames[0]}]} ${update_fnames[0]}| tail -n +2|cut -f ${update_ID_column[0]})"
    for i in $(seq 1 $((n_update-1)));do
	cmd=$cmd" <(${cats[${update_fnames[$i]}]} ${update_fnames[$i]}| tail -n +2|cut -f ${update_ID_column[$i]})"
    done
    cmn=$(eval $cmd | sort|uniq -c|sed 's/^  *//'|gawk -v n=$n_update '$1==n{print $0;}'| wc -l)
    echo "INFO: output $cmn common IDs between $n_update update files"|tee -a "$logfile"
    echo ""|tee -a "$logfile"

    echo -n "Merging update files ... "|tee -a "$logfile"

    # column classes of input columns are based on source update files
    for i in $(seq 0 $((n_update-1)));do
    	for c in $("${cats[${update_fnames[$i]}]}" "${update_fnames[$i]}"|head -n 1|cut --complement -f ${update_ID_column[$i]});do
    	    column_class[$c]=$i
    	done
    done
    column_class["$id_field"]="NA"

    if [[ $n_update -gt 1 ]];then
	join_cmd="join --header -t$'\t' -1 ${update_ID_column[0]} -2 ${update_ID_column[1]} <(cat <(${cats[${update_fnames[0]}]} ${update_fnames[0]}|head -n 1) <(${cats[${update_fnames[0]}]} ${update_fnames[0]}|tail -n +2|sort -k${update_ID_column[0]},${update_ID_column[0]})) <(cat <(${cats[${update_fnames[1]}]} ${update_fnames[1]}|head -n 1) <(${cats[${update_fnames[1]}]} ${update_fnames[1]}|tail -n +2|sort -k${update_ID_column[1]},${update_ID_column[1]}))"
	for i in $(seq 2 $((n_update-1)));do
    	    join_cmd=$join_cmd"|join --header -t$'\t' -1 1 -2 ${update_ID_column[$i]} - <(cat <(${cats[${update_fnames[$i]}]} ${update_fnames[$i]}|head -n 1) <(${cats[${update_fnames[$i]}]} ${update_fnames[$i]}|tail -n +2|sort -k${update_ID_column[$i]},${update_ID_column[$i]}))"
	done
	eval "$join_cmd > $tmpfile0"
    else
	cp "${update_fnames[$i]}" "$tmpfile0"
    fi
    
    head -n 1 "$tmpfile0" > "$tmpfile1"
    for c in $(head -n 1 "$tmpfile0");do
    	class_array+=(${column_class[$c]})
    done
    cline=$(join_by , "${class_array[@]}")
    echo $cline|tr ',' '\t' >> "$tmpfile1"
    tail -n +2 "$tmpfile0" >> "$tmpfile1"    
    rm -f "$tmpfile0"

    # # HEADER
    # command1="paste <(${cats[${update_fnames[0]}]} ${update_fnames[0]}|head -n 1)"
    # for i in $(seq 1 $((n_update-1)));do
    # 	command1=${command1}" <(${cats[${update_fnames[$i]}]} ${update_fnames[$i]}|head -n 1|cut --complement -f ${update_ID_column[$i]})"
    # done
    # # CLASS ROW
    # for c in $("${cats[${update_fnames[0]}]}" "${update_fnames[0]}"|head -n 1);do
    # 	class_array+=(${column_class[$c]})
    # done
    # for i in $(seq 1 $((n_update-1)));do
    # 	for c in $("${cats[${update_fnames[$i]}]}" "${update_fnames[$i]}"|head -n 1|cut --complement -f ${update_ID_column[$i]});do
    # 	    class_array+=(${column_class[$c]})
    # 	done
    # done
    # header2=$(join_by , "${class_array[@]}")
    # # BODY
    # command2="paste <(${cats[${update_fnames[0]}]} ${update_fnames[0]}|tail -n +2|sort -k${update_ID_column[0]},${update_ID_column[0]})"
    # for i in $(seq 1 $((n_update-1)));do
    # 	command2=${command2}" <(${cats[${update_fnames[$i]}]} ${update_fnames[$i]}|tail -n +2|sort -k${update_ID_column[$i]},${update_ID_column[$i]}|cut --complement -f ${update_ID_column[$i]})"
    # done
    # eval "cat <($command1) <(echo $header2|tr ',' '\t') <($command2) > $tmpfile1"
    
    echo "Done"|tee -a "$logfile"
    #----------------------------------------------------------------
    # input fields (except for ID column) and their classes, from the input file
    declare -a input_fields
    declare -a input_classes
    declare -A input_field2class
    for c in $("${cats[${input_fnames[0]}]}" "${input_fnames[0]}"|head -n 1|cut --complement -f ${input_ID_column[0]},${release_col},${created_col});do
	input_fields+=($c)
    done
    for c in $("${cats[${input_fnames[0]}]}" "${input_fnames[0]}"|head -n 2|tail -n 1|cut --complement -f ${input_ID_column[0]},${release_col},${created_col});do
	input_classes+=($c)
    done
    #echo "Input classes done"|tee -a "$logfile"
    for (( i=0; i<"${#input_fields[@]}"; i++ )); do input_field2class[${input_fields[$i]}]=${input_classes[$i]};done
    new_update_ID=$(getColNum "$tmpfile1" "$id_field" "cat")
    if [[ -z $new_update_ID ]];then
	echo "ERROR: no \"$id_field\" found in $tmpfile1"|tee -a "$logfile"
	exit 1
    fi    
    #echo "New update ID column: $new_update_ID"|tee -a "$logfile"
    # update fields (except for ID column) and their classes, from the merged update file
    declare -a update_fields
    declare -a update_classes
    declare -A update_field2class
    for c in $(head -n 1 "$tmpfile1"|cut --complement -f $new_update_ID);do update_fields+=($c);done
    #echo "Update fields done"|tee -a "$logfile"
    for c in $(head -n 2 "$tmpfile1"|tail -n 1|cut --complement -f $new_update_ID);do update_classes+=($c);done
    #echo "Update classes done"|tee -a "$logfile"
    for (( i=0; i<"${#update_fields[@]}"; i++ )); do update_field2class[${update_fields[$i]}]=${update_classes[$i]};done
    
    # input classes intersecting with update fields: all input fields from these will be excluded from input before merging
    declare -A intersecting_classes
    # new fields in update that will be added
    declare -A new_fields
    for c in $(head -n 1 "$tmpfile1"|cut --complement -f $new_update_ID);do
	# status=$(checkArray "$c" "${!input_field2class[@]}")
	# if [[ $status == "YES" ]];then
	status=${input_field2class[$c]}
	if [[ ! -z "$status" ]];then
	    x=${input_field2class[$c]}
	    intersecting_classes[$x]=1
	else
	    new_fields[$c]=1
	fi
    done
    #echo "Intersecting classes done"|tee -a "$logfile"
    
    # all column numbers from input intersecting_classes
    declare -a input_columns_to_exclude
    for (( i=0; i<"${#input_fields[@]}"; i++ )); do
	c=${input_fields[$i]}
	x=${input_field2class[$c]}
	# status=$(checkArray "$x" "${!intersecting_classes[@]}")
	# if [[ $status == "YES" ]];then
	status=${intersecting_classes[$x]}
	if [[ ! -z "$status" ]];then
	    n=$(getColNum ${input_fnames[0]} $c ${cats[${input_fnames[0]}]})
	    input_columns_to_exclude+=($n)
	fi
    done
    input_columns_to_exclude+=($created_col)
    input_columns_to_exclude+=($release_col)
    #echo "Columns to exclude done"|tee -a "$logfile"

    # input without excluded columns, without classes
    tmpfile2=$(mktemp -p "$out_dir" temp_2_input_XXXXXXXX)
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not create tmpfile2 "|tee -a "$logfile"
	exit 1
    fi
    
    echo -n "Removing columns from input ... "|tee -a "$logfile"
    #str=$(join_by , ${input_columns_to_exclude[*]})
    #echo "columns to exclude: $str"|tee -a "$logfile"
    "${cats[${input_fnames[0]}]}" "${input_fnames[0]}"|cut -f $(join_by , "${input_columns_to_exclude[@]}") --complement| gawk 'BEGIN{FS=OFS="\t";}NR!=2{print $0;}' > "$tmpfile2"
    echo "Done"|tee -a "$logfile"
    new_input_ID=$(getColNum "$tmpfile2" "$id_field" "cat")
    if [[ -z $new_input_ID ]];then
	echo "ERROR: no \"$id_field\" found in $tmpfile2"|tee -a "$logfile"
	exit 1
    fi    
    #echo "New input ID column: $new_input_ID"|tee -a "$logfile"

    common_IDs=$(join -1 1 -2 1 <(cut -f $new_update_ID "$tmpfile1"|tail -n +3|sort) <(cut -f $new_input_ID "$tmpfile2"|tail -n +2|sort)|wc -l)
    input_only_IDs=$(join -a 2 -e NULL -o 1.1,2.1 -1 1 -2 1 <(cut -f $new_update_ID "$tmpfile1"| tail -n +3|sort) <(cut -f $new_input_ID $tmpfile2|tail -n +2|sort)|grep NULL|wc -l)
    update_only_IDs=$(join -a 1 -e NULL -o 1.1,2.1 -1 1 -2 1 <(cut -f $new_update_ID "$tmpfile1"| tail -n +3|sort) <(cut -f $new_input_ID $tmpfile2|tail -n +2|sort)|grep NULL|wc -l)
    echo ""|tee -a "$logfile"
    echo "INFO: common IDs between input and update: $common_IDs"|tee -a "$logfile"
    echo "INFO: IDs only in input: $input_only_IDs (not included in output)"|tee -a "$logfile"
    str=""
    if [[ $mode == "inner" ]];then
	str=" (not included in output, mode=$mode)"
    else
	str=" (included in output, mode=$mode)"
    fi
    echo "INFO: IDs only in update: ${update_only_IDs}$str"|tee -a "$logfile"
    echo ""|tee -a "$logfile"

    fmtstr="1.""${new_input_ID}"",2.""${new_update_ID}"
    input_ncol=$(head -n 1 "$tmpfile2"|tr '\t' '\n'|wc -l)
    update_ncol=$(head -n 1 "$tmpfile1"|tr '\t' '\n'|wc -l)
    for (( i=1; i<=$input_ncol; i++ )); do
	if [[ $i -ne $new_input_ID ]];then
	    fmtstr=$fmtstr",1.$i"
	fi
    done
    for (( i=1; i<=$update_ncol; i++ )); do
	if [[ $i -ne $new_update_ID ]];then
	    fmtstr=$fmtstr",2.$i"
	fi
    done
    
    # join input and merged update, without classes
    tmpfile3=$(mktemp -p "$out_dir" temp_3_joined_XXXXXXXX)
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not create tmpfile3 "|tee -a "$logfile"
	exit 1
    fi
    
    echo -n "Joining ... "|tee -a "$logfile"
    join --header -1 $new_input_ID -2 $new_update_ID -a 1 -a 2 -t $'\t' -e NA -o $fmtstr "$tmpfile2" <(gawk 'BEGIN{FS=OFS="\t";}NR!=2{print $0;}' "$tmpfile1")|gawk -v m=$mode 'BEGIN{FS=OFS="\t";}{if (NR==1){print $0;}else{if ($2!="NA"){if (m=="inner"){if ($1!="NA"){print $0;}} if(m=="right"){if($1=="NA"){$1=$2;}print $0;} } }}'| cut -f 2 --complement > "$tmpfile3"
    echo "Done"|tee -a "$logfile"
    joined_ID=$(getColNum $tmpfile3 "$id_field" "cat")
    if [[ -z $joined_ID ]];then
	echo "ERROR: no \"$id_field\" found in $tmpfile3"|tee -a "$logfile"
	exit 1
    fi    
    #echo "Joined file input ID column: $joined_ID"|tee -a "$logfile"

    # max update class
    maxc=0
    for c in "${update_field2class[@]}";do
	if [[ $c -gt $maxc ]];then
	    maxc=$c
	fi
    done
    declare -A new_classes
    declare -A tmp_ar
    cur_add=1
    for c in $(head -n 1 "$tmpfile3"|cut -f $joined_ID --complement);do
	# s=$(checkArray "$c" "${!update_field2class[@]}")
	# if [[ $s == "YES" ]];then
	s=${update_field2class[$c]}
	if [[ ! -z "$s" ]];then
	    new_classes[$c]=${update_field2class[$c]}
	else
	    # t=$(checkArray "$c" "${!input_field2class[@]}")
	    # if [[ $t == "YES" ]];then
	    t=${input_field2class[$c]}
	    if [[ ! -z "$t" ]];then
		d=${input_field2class[$c]}
		# u=$(checkArray "$d" "${!tmp_ar[@]}")
		# if [[ $u == "YES" ]];then
		u=${tmp_ar[$d]}
		if [[ ! -z "$u" ]];then
		    new_classes[$c]=${tmp_ar[$d]}
		else
		    tmp_ar[$d]=$((maxc+cur_add))
		    new_classes[$c]=${tmp_ar[$d]}
		    cur_add=$((cur_add+1))
		fi
	    else
		echo "ERROR: could not find class for $c"|tee -a "$logfile"
		exit 1
	    fi
	fi
    done
    #echo "New classes done"|tee -a "$logfile"

    # FINAL OUTPUT
    echo -n "Creating output file ... "|tee -a "$logfile"
    command1="paste <(head -n 1 $tmpfile3) <(echo RELEASE CREATED|tr ' ' '\t')"

    declare -a temp
    for c in $(head -n 1 "$tmpfile3");do
	if [[ $c == "$id_field" ]];then
	    temp+=("NA")
	else
	    temp+=(${new_classes[$c]})
	fi
    done
    temp+=("NA")
    temp+=("NA")
    header2=$(join_by , "${temp[@]}")

    command2="paste <(tail -n +2 $tmpfile3)"
    x=$(tail -n +2 "$tmpfile3"| wc -l)
    command2=${command2}" <(yes $release $datestr| tr ' ' '\t'| head -n $x)"

    eval "cat <($command1) <(echo $header2|tr ',' '\t') <($command2) | gzip - -c > ${outfile}"
    echo "Done"|tee -a "$logfile"

    rm -f "$tmpfile1"
    rm -f "$tmpfile2"
    rm -f "$tmpfile3"
    
    echo "Done"|tee -a "$logfile"
    date "+%F %H-%M-%S" |tee -a "$logfile"
fi

exit 0

