#!/usr/bin/env bash

scriptdir=$(dirname $(readlink -f $0))
upperdir=$(dirname $scriptdir)
source "${upperdir}/functions.sh"
scriptname="${scriptdir}/selectCases2.pl"

nsamples=$1
out_table=$2
out_codes=$3
out_cases=$4

declare -a codes
declare -a rules
declare -a expressions
declare -a temp
declare -a cases

max_inst=10
min_inst=2
max_codes=5
nrules=100

# ------------------------------------------------------------------------------------------------------------------

for c in {A..B};do
    for i in $(seq -w 0 999);do
	codes+=($c$i)
    done
done
codes_str=$(join_by , "${codes[@]}")

while read e;do
    expressions+=("$e")
done < <(perl -le '%a=("X"=>1,"X and X"=>1,"X or X"=>1);for ($i=0;$i<3;$i++){$k=(keys %a)[rand keys %a];if ($k ne "X"){$c=()=$k=~/X/g;for ($j=1;$j<=$c;$j++){$x=(keys %a)[rand keys %a];next if ($x eq "X");$z=$j;$k=~s/(X)/--$z==0?"(".$x.")":$1/ge;$a{$k}=1;}}}foreach $x (keys %a){print $x;}')
exp_str=$(join_by , "${expressions[@]}")

:> "$out_codes"
for i in $(seq 1 $nrules);do
    e=$(echo $exp_str|tr ',' '\n'|shuf -n 1)
    temp=()
    x=$(echo "$e"|sed 's/[^X]//g'|awk '{print length;}')
    while read c;do
	temp+=($c)
    done < <(echo $codes_str|tr ',' '\n'|shuf -n $x)
    str=$(join_by , "${temp[@]}")
    s=$(perl -se 'foreach $x (split(/,/,$b,-1)){$a=~s/X/$x/;}print $a."\n";' -- -a="$e" -b="$str")
    echo $s >> "$out_codes"
    rules+=("$s")
done

echo "eid" "ins_index" "arr_index" "diag_icd10" "diag_icd9" | tr ' ' '\t' > "$out_table"
for id in $(seq -w 1 $nsamples);do
    inst=$((RANDOM%max_inst))
    inst=$((inst+min_inst))
    for i in $(seq 1 $inst);do
	temp=()
	x=$((RANDOM%max_codes))
	x=$((x+1))
	j=1
	while read c;do
	    temp+=($c)
	    echo $id $i $j $c "ICD9" | tr ' ' '\t' >> "$out_table"
	    j=$((j+1))
	done < <(echo $codes_str|tr ',' '\n'|shuf -n $x)
      	str=$(join_by , "${temp[@]}")	
	echo $id $i $str| tr ' ' '\t'
    done
done | "$scriptname" "$out_codes" 2|grep -v False > "$out_cases"

exit 0
