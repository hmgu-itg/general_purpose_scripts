#!/usr/bin/env bash

scriptname=$0
args=("$@")

declare -i rep=$1
declare -i rows=$2
declare -i cols=$3
declare -i n=$4
outdir=$5

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
ginput="${scriptdir}/generate_input.sh"
gupdate="${scriptdir}/generate_update.sh"
mergesh="${scriptdir}/ukbb_merge.sh"
mergepy="${scriptdir}/ukbb_merge.py"

tmpdir=$(mktemp -d -p "$outdir" tempdir_test_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary directory in $outdir"
    exit 1
fi

"$ginput" $rows $cols $n "$tmpdir"/input
echo ""

for (( k=1; k <= $rep; k++ ));do
    echo "------------------------------------ ITERATION $k ------------------------------------"
    echo ""
    opts=""
    if [[ $k -eq 1 ]];then
	for (( i=1; i <= $n; i++ ));do
	    opts="$opts -i $tmpdir/input_$i.txt.gz" 
	done
    else
	zcat "$infile"|cut -f 1|tail -n +3|shuf -n 10 > "${tmpdir}"/exclude
	opts="-i $infile -x ${tmpdir}/exclude"
	for (( i=1; i <= $n; i++ ));do
	    opts="$opts -u ${tmpdir}/update_r${k}_${i}.txt.gz" 
	done
    fi

    echo "OPTS: $opts"
    echo ""

    eval "$mergesh" "$opts" -b merge_sh -o "$tmpdir" 1>/dev/null
    eval "$mergepy" "$opts" -o "$tmpdir"/merge_py_r$k
    zdiff "$tmpdir"/merge_sh_r$k.txt.gz "$tmpdir"/merge_py_r$k.txt.gz &>/dev/null
    if [[ $? -ne 0 ]];then
	echo "ERROR: iteration $k, files $tmpdir/merge_sh_r$k.txt.gz $tmpdir/merge_py_r$k.txt.gz differ"
	exit 1
    else
	echo ""
	echo "INFO: iteration $k test OK"
	echo ""
    fi

    infile="$tmpdir"/merge_sh_r$k.txt.gz
    echo "UPDATING"
    echo ""
    "$gupdate" $infile "$tmpdir"/update >> "$tmpdir"/update.log
    echo ""
done

exit 0

