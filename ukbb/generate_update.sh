#!/usr/bin/env bash

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

infile=$1
prefix=$2

catcmd=$(getCatCmd $infile)
checkFields $infile $catcmd

declare -a classes

ns=$($catcmd $infile|tail -n +3| wc -l)
ncols=$($catcmd $infile|head -n 1|tr '\t' '\n'|wc -l)
ncols=$((ncols-3))

max_mod=$((ns/100))

if [[ $max_mod -eq 0 ]];then max_mod=1;fi

del_col_lim=3
add_col_lim=3

idcol=$(getColNum $infile "f.eid" $catcmd)
relcol=$(getColNum $infile "RELEASE" $catcmd)
crcol=$(getColNum $infile "CREATED" $catcmd)
if [[ -z "$idcol" ]];then
    echo "ERROR: ID column not found in $infile"
    exit 1
fi
if [[ -z "$relcol" ]];then
    echo "ERROR: RELEASE column not found in $infile"
    exit 1
fi
if [[ -z "$crcol" ]];then
    echo "ERROR: CREATED column not found in $infile"
    exit 1
fi

echo "ID column: $idcol"
echo "RELEASE column: $relcol"
echo "CREATED column: $crcol"

release=$($catcmd $infile|cut -f $relcol|head -n 3|tail -n 1)
release=$((release+1))
echo "New release: $release"
echo ""

read -r -a classes <<<$($catcmd $infile|head -n 2|tail -n 1|cut --complement -f ${idcol},${relcol},${crcol}|tr '\t' '\n'|sort -n|uniq|tr '\n' ' '|sed 's/ $//')

declare -a sar
i=1
newcol_index_start=1
newcol_index_end=""
for c in "${classes[@]}";do
    echo "Generating file $i/${#classes[@]}"
    fname=${prefix}_r${release}_${i}".txt.gz"
    echo "Output file: $fname"

    to_mod=$((RANDOM%max_mod))
    rem=$((ns-to_mod))
    del=$((RANDOM%2))

    if [[ $to_mod -eq 0 ]];then
	id_cmd='$catcmd $infile|cut -f 1|tail -n +3'
	new_total=$ns
	echo "Keeping same IDs"
    else
	if [[ $del -eq 1 ]];then
	    echo "Deleting $to_mod row(s)"
	    tmp=$($catcmd $infile|cut -f 1|tail -n +3|shuf|head -n $rem|sort|tr '\n' ',')
	    id_cmd='echo -n '$tmp'|sed "s/,/\n/g"'
	    new_total=$((ns-to_mod))
	else
	    echo "Adding $to_mod row(s)"
	    id_cmd='cat <($catcmd $infile|cut -f 1|tail -n +3) <(for i in $(seq -w 1 '$to_mod');do echo ADD$i;done)|sort'
	    new_total=$((ns+to_mod))
	fi
    fi
    
    read -r -a sar <<<$($catcmd $infile|head -n 2|tail -n 1|tr '\t' '\n'|cat -n|sed 's/^  *//'|sed 's/\t/ /g'|gawk -v n=$c '$2==n{print $1;}'|tr '\n' ' '|sed 's/ $//')
    str=$(join_by , "${sar[@]}")
#    echo $str
    n="${#sar[@]}"
    fmt="1.1"
    for (( j=0; j<$n; j++ ));do
	fmt=$fmt",2.$((j+2))"
    done
#    echo $fmt
    to_del=$((RANDOM%del_col_lim))
    echo "Deleting $to_del column(s)"
    to_del_str=""
    if [[ $to_del -gt 0 ]];then
	to_del_str=$(seq 2 $((n+1))|shuf|head -n $to_del|sort -n|tr '\n' ','|sed 's/,$//')
#	echo "to_del_str=${to_del_str}"
    fi
    to_add=$((RANDOM%add_col_lim))
    echo "Adding $to_add column(s)"
    echo ""

    cmd=""
    if [[ -z "${to_del_str}" ]];then
	cmd="join -t $'\t' --header -1 1 -2 1 -a 1 -a 2 -e NULL -o $fmt <(cat <(echo f.eid) <(eval $id_cmd)) <(cat <($catcmd $infile|head -n 1|cut -f $idcol,$str) <($catcmd $infile|tail -n +3|cut -f $idcol,$str|sort -k1,1)) | grep -v ^NULL"
    else
	cmd="join -t $'\t' --header -1 1 -2 1 -a 1 -a 2 -e NULL -o $fmt <(cat <(echo f.eid) <(eval $id_cmd)) <(cat <($catcmd $infile|head -n 1|cut -f $idcol,$str) <($catcmd $infile|tail -n +3|cut -f $idcol,$str|sort -k1,1)) | grep -v ^NULL|cut --complement -f ${to_del_str}"
    fi

    if [[ $to_add -gt 0 ]];then
	newcol_index_end=$((newcol_index_start+to_add-1))
	cmd2="cat <(seq -w -s' ' $newcol_index_start $newcol_index_end|sed -e 's/^/NEWCOL_${release}_/g' -e 's/ / NEWCOL_${release}_/g' -e 's/ /\t/g') <(yes $(seq -s ' ' 1 $to_add|sed 's/[0-9][0-9]*/NEWVAL/g')|head -n $new_total|tr ' ' '\t')"
	cmd="paste <($cmd) <($cmd2)"
	newcol_index_start=$((newcol_index_end+1))
    fi
    
    eval "$cmd | gzip - -c > $fname"
    
    i=$((i+1))
done
