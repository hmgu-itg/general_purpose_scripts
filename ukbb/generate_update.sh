#!/usr/bin/env bash

function getColNum () {
    local fname=$1
    local colname=$2
    echo $(fgrep -w $colname  <(head -n 1 $fname | tr '\t' '\n'| cat -n | sed 's/^  *//') | cut -f 1)
}

function join_by { local IFS="$1"; shift; echo "$*"; }

infile=$1
prefix=$2

declare -a classes

ns=$(tail -n +3 $infile| wc -l)
ncols=$(head -n 1 $infile|tr '\t' '\n'|wc -l)
ncols=$((ncols-3))

max_mod=$((ns/100))

if [[ $max_mod -eq 0 ]];then max_mod=1;fi
to_mod=$((RANDOM%max_mod))
rem=$((ns-to_mod))
del=$((RANDOM%2))

echo "to_mod=$to_mod" 1>&2

if [[ $to_mod -eq 0 ]];then
    id_cmd='cut -f 1 $infile|tail -n +3'
    new_total=$ns
else
    if [[ $del -eq 1 ]];then
	echo "Deleting $to_mod row(s)" 1>&2
	tmp=$(cut -f 1 $infile|tail -n +3|shuf|head -n $rem|sort|tr '\n' ',')
	id_cmd='echo -n '$tmp'|sed "s/,/\n/g"'
	new_total=$((ns-to_mod))
    else
	echo "Adding $to_mod row(s)" 1>&2
	id_cmd='cat <(cut -f 1 '$infile'|tail -n +3) <(for i in $(seq -w 1 '$to_mod');do echo ADD$i;done)|sort'
	new_total=$((ns+to_mod))
    fi
fi
echo ""

del_col_lim=3
add_col_lim=3

idcol=$(getColNum $infile "f.eid")
relcol=$(getColNum $infile "RELEASE")
crcol=$(getColNum $infile "CREATED")

read -r -a classes <<<$(head -n 2 $infile|tail -n 1|cut --complement -f ${idcol},${relcol},${crcol}|tr '\t' '\n'|sort -n|uniq|tr '\n' ' '|sed 's/ $//')

declare -a sar
i=1
for c in ${classes[@]};do
    echo "Generating file $i/${#classes[@]}" 1>&2
    fname=${prefix}_${i}".txt"
    echo "Output file: $fname" 1>&2
    
    read -r -a sar <<<$(head -n 2 $infile|tail -n 1|tr '\t' '\n'|cat -n|sed 's/^  *//'|sed 's/\t/ /g'|awk -v n=$c '$2==n{print $1;}'|tr '\n' ' '|sed 's/ $//')
    str=$(join_by , ${sar[*]})
    echo $str 1>&2
    n=${#sar[@]}
    fmt="1.1"
    for (( j=0; j<$n; j++ ));do
	fmt=$fmt",2.$((j+2))"
    done
    echo $fmt 1>&2
    to_del=$((RANDOM%del_col_lim))
    echo "Deleting $to_del column(s)"
    to_del_str=""
    if [[ $to_del -gt 0 ]];then
	to_del_str=$(seq 2 $((n+1))|shuf|head -n $to_del|sort -n|tr '\n' ','|sed 's/,$//')
	echo "to_del_str=${to_del_str}"
    fi
    to_add=$((RANDOM%add_col_lim))
    echo "Adding $to_add column(s)"
    echo ""

    cmd=""
    if [[ -z "${to_del_str}" ]];then
	cmd="join -t $'\t' --header -1 1 -2 1 -a 1 -a 2 -e NULL -o $fmt <(cat <(echo f.eid) <(eval $id_cmd)) <(cat <(head -n 1 $infile|cut -f $idcol,$str) <(tail -n +3 $infile|cut -f $idcol,$str|sort -k1,1)) | grep -v ^NULL"
    else
	cmd="join -t $'\t' --header -1 1 -2 1 -a 1 -a 2 -e NULL -o $fmt <(cat <(echo f.eid) <(eval $id_cmd)) <(cat <(head -n 1 $infile|cut -f $idcol,$str) <(tail -n +3 $infile|cut -f $idcol,$str|sort -k1,1)) | grep -v ^NULL|cut --complement -f ${to_del_str}"
    fi

    if [[ $to_add -gt 0 ]];then
	cmd2="cat <(seq -w -s' ' 1 $to_add|sed -e 's/^/NEWCOL/g' -e 's/ / NEWCOL/g' -e 's/ /\t/g') <(yes $(seq -s ' ' 1 $to_add|sed 's/[0-9][0-9]*/NEWVAL/g')|head -n $new_total|tr ' ' '\t')"
	cmd="paste <($cmd) <($cmd2)"
    fi
    
    eval "$cmd > $fname"
    
    i=$((i+1))
done
