#!/usr/bin/env bash

# ALL FILES ARE TAB SEPARATED
# for all fields (except f.eid), report files where they occur

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

function usage () {
    echo ""
    echo "Report fields overlap"
    echo ""
    echo "Usage: fields_overlap.sh -f <ID field name; default: \"f.eid\">"
    echo "                         -i input1.tab"
    echo "                         -i input2.tab"
    echo "           ...                        "
    echo "                         -i inputN.tab"
    echo ""
    exit 0
}

declare -a input_fnames
declare -a bnames
declare -a input_ID_column
declare -A cats

if [[ $# -eq 0 ]];then
    usage
fi

id_field="f.eid"
while getopts "hi:f:" opt; do
    case $opt in
        i)input_fnames+=($OPTARG);;
        f)id_field=($OPTARG);;
        h)usage;;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

n_input="${#input_fnames[@]}"

if [[ $n_input -eq 0 ]];then
    echo "ERROR: no input files specified" 1>&2
    exit 1
fi

# get zcat/cat command for each input file
for i in $(seq 0 $((n_input-1)));do
    cats["${input_fnames[$i]}"]=$(getCatCmd "${input_fnames[$i]}")
done

#----------------------------------------------------

# get ID column number for every input file
for i in $(seq 0 $((n_input-1)));do
    x=$(getColNum "${input_fnames[$i]}" "$id_field" ${cats[${input_fnames[$i]}]})
    exitIfEmpty "$x" "ERROR: no \"$id_field\" found in ${input_fnames[$i]}"
    input_ID_column+=($x)
    bnames+=($(basename ${input_fnames[$i]}))
done

#----------------------------------------------------

declare -A header_arr
for i in $(seq 0 $((n_input-1)));do
    >&2 echo "INFO: checking file header $i"
    b="${bnames[$i]}"
    while read cname;do
	if [[ -v "header_arr[$cname]" ]];then
	    header_arr[$cname]="${header_arr[$cname]}"",$b"
	else
	    header_arr[$cname]="$b"
	fi
    done < <(${cats["${input_fnames[$i]}"]} "${input_fnames[$i]}" | head -n 1 | cut --complement -f "${input_ID_column[$i]}"| tr '\t' '\n')
done

for cname in "${!header_arr[@]}";do
    echo $cname ${header_arr[$cname]} | tr ' ' '\t'
done
    
exit 0

