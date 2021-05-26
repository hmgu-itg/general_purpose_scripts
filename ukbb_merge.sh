#!/usr/bin/env bash

# ALL INPUT FILES ARE TSV

function getColNum () {
    local fname=$1
    local colname=$2
    echo $(fgrep -w $colname  <(head -n 1 $fname | tr '\t' '\n'| cat -n | sed 's/^  *//') | cut -f 1)
}

function usage () {
    echo "Merging script"
}

# default ID column
id_field="f.eid"

declare -a input_fnames
declare -a update_filenames
declare -a input_ID_column
declare -a update_ID_column
declare -a input_nrows
declare -a update_nrows

while getopts "i:u:f:" opt; do
    case $opt in
        i)input_fnames+=($OPTARG);;
        u)update_fnames+=($OPTARG);;
        f)id_field=($OPTARG);;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

n_input=${#input_fnames[@]}
n_update=${#update_fnames[@]}

if [[ $n_input -eq 0 ]];then
    echo "ERROR: no input files specified" 1>&2
    exit 1;
fi

if [[ $n_input -lt 2 && $n_update -eq 0 ]];then
    echo "ERROR: no update files specified, only one input file specified (need at least two input files to merge)" 1>&2
    exit 1;
fi

if [[ $n_update -gt 0 && $n_input -gt 1 ]];then
    echo "WARN: $n_update update files and $n_input input files specified; only the first input file (${input_fnames[0]}) will be processed" 1>&2
fi

for i in $(seq 0 $((n_input-1)));do
    x=$(getColNum ${input_fnames[$i]} $id_field)
    input_ID_column+=($x)
    input_nrows+=$(cat ${input_fnames[$i]} | wc -l)
done

for i in $(seq 0 $((n_update-1)));do
    x=$(getColNum ${update_fnames[$i]} $id_field)
    update_ID_column+=($x)
    update_nrows+=$(cat ${update_fnames[$i]} | wc -l)
done

# output info
for i in $(seq 0 $((n_input-1)));do
    echo "INFO: input file: ${input_fnames[$i]}" 1>&2
    echo "INFO: ID column index: ${input_ID_column[$i]}" 1>&2
    echo "INFO: rows: ${input_nrows[$i]}" 1>&2
    echo "" 1>&2
done

for i in $(seq 0 $((n_update-1)));do
    echo "INFO: update file: ${update_fnames[$i]}" 1>&2
    echo "INFO: ID column index: ${update_ID_column[$i]}" 1>&2
    echo "INFO: rows: ${update_nrows[$i]}" 1>&2
    echo "" 1>&2
done

    

# nrows=${input_nrows[0]}
# for i in $(seq 0 $((n_input-1)));do
#     x=${input_nrows[$i]}
#     if [[ $x -ne $nrows ]];then
# 	echo "ERROR: files ${input_fnames[0]}  and ${input_fnames[$i]} have different number of rows" 1>&2
# 	exit 1
#     fi
# done

# for i in $(seq 0 $((n_update-1)));do
#     x=${update_nrows[$i]}
#     if [[ $x -ne $nrows ]];then
# 	echo "ERROR: files ${input_fnames[0]}  and ${update_fnames[$i]} have different number of rows" 1>&2
# 	exit 1
#     fi
# done

#
# check if all files have the same sample IDs
#
echo "Checking if input/update files have same sample IDs ... " 1>&2
for i in $(seq 0 $((n_input-1)));do
    x=$(cat <(cut -f ${input_ID_column[0]} ${input_fnames[0]}) <(cut -f ${input_ID_column[$i]} ${input_fnames[$i]})|sort|uniq -u|wc -l)
    if [[ $x -ne 0 ]];then
	echo "ERROR: files ${input_fnames[0]} and ${input_fnames[$i]} have different sets of IDs" 1>&2
	exit 1
    fi
done

for i in $(seq 0 $((n_update-1)));do
    x=$(cat <(cut -f ${input_ID_column[0]} ${input_fnames[0]}) <(cut -f ${update_ID_column[$i]} ${update_fnames[$i]})|sort|uniq -u|wc -l)
    if [[ $x -ne 0 ]];then
	echo "ERROR: files ${input_fnames[0]} and ${update_fnames[$i]} have different sets of IDs" 1>&2
	exit 1
    fi
done
echo "OK" 1>&2
echo "" 1>&2

#
# check if column names (excepth for ID field) in input files are disjoint
#
echo "Checking if column names in input/update files are disjoint ... " 1>&2
if [[ $n_input -gt 1 ]];then
    for i in $(seq 0 $((n_input-1)));do
	for j in $(seq $((i+1)) $((n_input-1)));do
	    x=$(cat <(head -n 1 ${input_fnames[$i]}|cut --complement -f ${input_ID_column[$i]}) <(head -n 1 ${input_fnames[$j]}|cut --complement -f ${input_ID_column[$j]})|sort|uniq -d|wc -l)
	    if [[ $x -ne 0 ]];then
		echo "ERROR: input files ${input_fnames[$i]} and ${input_fnames[$j]} have columns in common" 1>&2
		exit 1
	    fi	    
	done
    done
fi

#
# check if column names (excepth for ID field) in update files are disjoint
#
if [[ $n_update -gt 1 ]];then
    for i in $(seq 0 $((n_update-1)));do
	for j in $(seq $((i+1)) $((n_update-1)));do
	    x=$(cat <(head -n 1 ${update_fnames[$i]}|cut --complement -f ${update_ID_column[$i]}) <(head -n 1 ${update_fnames[$j]}|cut --complement -f ${update_ID_column[$j]})|sort|uniq -d|wc -l)
	    if [[ $x -ne 0 ]];then
		echo "ERROR: update files ${update_fnames[$i]} and ${update_fnames[$j]} have columns in common" 1>&2
		exit 1
	    fi	    
	done
    done
fi
echo "OK" 1>&2
echo "" 1>&2

if [[ $n_update -eq 0 ]];then # just merging input files
    echo "Merging input files ... " 1>&2
    command1="paste <(head -n 1 ${input_fnames[0]})"
    for i in $(seq 1 $((n_input-1)));do
	command1=$command1" <(head -n 1 ${input_fnames[$i]}|cut --complement -f ${input_ID_column[$i]})"
    done

    command2="paste <(tail -n +2 ${input_fnames[0]}|sort -k${input_ID_column[0]},${input_ID_column[0]})"
    for i in $(seq 1 $((n_input-1)));do
	command2=$command2" <(tail -n +2 ${input_fnames[$i]}|sort -k${input_ID_column[$i]},${input_ID_column[$i]}|cut --complement -f ${input_ID_column[$i]})"
    done

    eval "cat <($command1) <($command2)"
    echo "Done" 1>&2
    echo "" 1>&2
else # updating the first input file using update files
    exit 0
    echo "Merging input files ... " 1>&2

    echo "Done" 1>&2
    echo "" 1>&2

    
fi

# cat <(paste <(head -n 1 $fname1) <(head -n 1 $fname2| cut -f 2-) <(head -n 1 $fname3| cut -f 2-)) <(paste -d ',' <(tail -n +2 $fname1|sort -t ',' -k1,1) <(tail -n +2 $fname2| sort -t ',' -k1,1|cut -d ',' -f 2-) <(tail -n +2 $fname3|sort -t ',' -k1,1| cut -d ',' -f 2-))

exit 0

