#!/usr/bin/env bash

# ALL INPUT FILES ARE TAB SEPARATED

bold=$(tput bold)
underlined=$(tput smul)
normal=$(tput sgr0)

function checkFields {
    local fname=$1
    echo -n "Checking # fields in $fname ... "
    local x=$(awk 'BEGIN{FS="\t";}{print NF;}' $fname| sort|uniq| wc -l)
    if [[ $x -eq 1 ]];then
	echo "OK"
    else
	echo "ERROR: $fname contains rows with different number of fields" 1>&2
	exit 1
    fi
}

function getColNum () {
    local fname=$1
    local colname=$2
    echo $(fgrep -w $colname  <(head -n 1 $fname | tr '\t' '\n'| cat -n | sed 's/^  *//') | cut -f 1)
}

function usage () {
    echo ""
    echo "Merging script for inpatient UKBB data"
    echo ""
    echo "Usage: hesin_merge.sh -i <main HESIN table> -d <HESIN DIAG table> -p <HESIN OPER table>"
    echo ""
    echo "All input/update files are tab-separated."
    echo ""
    exit 0
}

while getopts "hi:d:p:" opt; do
    case $opt in
        f)main_fname=($OPTARG);;
        f)diag_fname=($OPTARG);;
        f)oper_fname=($OPTARG);;
        h)usage;;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

checkFields $main_fname
checkFields $diag_fname
checkFields $oper_fname

main_cols=$(head -n 1 $main_fname|tr '\t' '\n'|wc -l)
diag_cols=$(head -n 1 $diag_fname|tr '\t' '\n'|wc -l)
oper_cols=$(head -n 1 $oper_fname|tr '\t' '\n'|wc -l)

main_str="1.1"
for i in $(seq 2 $((main_cols-1)));do
    main_str=$main_str",1.$i"
done

diag_str="2.2"
for i in $(seq 3 $((diag_cols-1)));do
    diag_str=$diag_str",2.$i"
done


# IDs are first 2 fields

join --header -1 1 -2 1 -e NA -a 1 -a 2 -o ${main_str},${diag_str} <(cat <(head -n 1 $main_fname| sed 's/\t/\./') <(tail -n +2 $main_fname| sed 's/\t/\./'| sort -k1,1)) <(cat <(head -n 1 $diag_fname| sed 's/\t/\./') <(tail -n +2 $diag_fname| sed 's/\t/\./'| sort -k1,1)) | tr '\.' ',' | tr ' ' ','
