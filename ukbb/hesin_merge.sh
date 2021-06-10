#!/usr/bin/env bash

# ALL INPUT FILES ARE TAB SEPARATED

bold=$(tput bold)
underlined=$(tput smul)
normal=$(tput sgr0)

source "functions.sh"

function usage () {
    echo ""
    echo "Merging script for inpatient UKBB data"
    echo ""
    echo "Usage: hesin_merge.sh -i <main HESIN table> -d <HESIN DIAG table> -p <HESIN OPER table> -r <release>"
    echo ""
    echo "All input/update files are tab-separated."
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

while getopts "hi:d:p:r:" opt; do
    case $opt in
        i)main_fname=($OPTARG);;
        d)diag_fname=($OPTARG);;
        p)oper_fname=($OPTARG);;
        r)release=($OPTARG);;
        h)usage;;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

echo "MAIN TABLE: $main_fname" 1>&2
echo "DIAG TABLE: $diag_fname" 1>&2
echo "OPER TABLE: $oper_fname" 1>&2
echo "" 1>&2

checkFields $main_fname
checkFields $diag_fname
checkFields $oper_fname
echo "" 1>&2

main_cols=$(head -n 1 $main_fname|tr '\t' '\n'|wc -l)
diag_cols=$(head -n 1 $diag_fname|tr '\t' '\n'|wc -l)
oper_cols=$(head -n 1 $oper_fname|tr '\t' '\n'|wc -l)

echo "COLUMNS IN MAIN TABLE: $main_cols" 1>&2
echo "COLUMNS IN DIAG TABLE: $diag_cols" 1>&2
echo "COLUMNS IN OPER TABLE: $oper_cols" 1>&2
echo "" 1>&2

main_str="1.1"
for i in $(seq 2 $((main_cols-1)));do
    main_str=$main_str",1.$i"
done

diag_str="2.2"
for i in $(seq 3 $((diag_cols-1)));do
    diag_str=$diag_str",2.$i"
done

join_str="1.1"
for i in $(seq 2 $((main_cols+diag_cols-3)));do
    join_str=$join_str",1.$i"
done

oper_str="2.2"
for i in $(seq 3 $((oper_cols-1)));do
    oper_str=$oper_str",2.$i"
done

echo "MAIN STR: $main_str" 1>&2
echo "DIAG STR: $diag_str" 1>&2
echo "JOIN STR: $join_str" 1>&2
echo "OPER STR: $oper_str" 1>&2

datestr=$(date +%F)

# IDs are the first 2 fields

join --header -t $'\t' -1 1 -2 1 -e NA -a 1 -a 2 -o ${main_str},${diag_str} <(cat <(head -n 1 $main_fname| sed 's/\t/\./'|awk 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i ~ /^ *$/){$i="NA";}}print $0;}') <(tail -n +2 $main_fname| sed 's/\t/\./'| sort -k1,1|awk 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i ~ /^ *$/){$i="NA";}}print $0;}')) <(cat <(head -n 1 $diag_fname|sed 's/level/diag_level/'|sed 's/arr_index/diag_arr_index/'| sed 's/\t/\./'|awk 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i ~ /^ *$/){$i="NA";}}print $0;}') <(tail -n +2 $diag_fname| sed 's/\t/\./'| sort -k1,1|awk 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i ~ /^ *$/){$i="NA";}}print $0;}'))| join --header -t $'\t' -1 1 -2 1 -e NA -a 1 -a 2 -o ${join_str},${oper_str} - <(cat <(head -n 1 $oper_fname|sed 's/level/oper_level/'|sed 's/arr_index/oper_arr_index/'|sed 's/\t/\./'|awk 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i ~ /^ *$/){$i="NA";}}print $0;}') <(tail -n +2 $oper_fname| sed 's/\t/\./'| sort -k1,1|awk 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i ~ /^ *$/){$i="NA";}}print $0;}'))|sed 's/\./\t/' | awk -v d=$datestr -v r=$release 'BEGIN{FS=OFS="\t";}{if (NR==1){print $0,"CREATED","RELEASE";}else{print $0,d,r;}}'
