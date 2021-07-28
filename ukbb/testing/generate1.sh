#!/usr/bin/env bash

scriptdir=$(dirname $(readlink -f $0))
upperdir=$(dirname $scriptdir)
source "${upperdir}/functions.sh"

outname1=$1
outname2=$2

declare -ir max_arrays=10
declare -ir max_visits=10
declare -r p1=0.5
declare -a opcodes
declare -a icd10codes
declare -a icd9codes
declare -ir out_opcodes=50
declare -r nsamples=1000

# perl -e 'sub F{my $l=shift;my @a=("AND","OR");my $p=rand(2);return "x ".$a[$p]." x" if $l==5;my $x1=F($l+1);my $x2=F($l+1);return "(".$x1.") ".$a[$p]." (".$x2.")";} print F(1)."\n";'

declare -a expressions=("X And (X Or X Or X Or (X AND X) Or X)" "X And (X Or X Or X)")

# ------------------------------------------------------------------------------------------------------------------

for c in {A..B};do
    for i in $(seq -w 0 999);do
	opcodes+=($c$i)
    done
done
opcodes_str=$(join_by , "${opcodes[@]}")

for c in {K..N};do
    for i in $(seq -w 0 99);do
	icd10codes+=($c$i)
    done
done
icd10codes_str=$(join_by , "${icd10codes[@]}")

for c in {Q..T};do
    for i in $(seq -w 0 99);do
	icd9codes+=($c$i)
    done
done
icd9codes_str=$(join_by , "${icd9codes[@]}")

op_sz="${#opcodes[@]}"
sz2="${#expressions[@]}"

for (( i=0; i<$out_opcodes; i++ ));do
    x=$((RANDOM%10))
    if [[ $x -lt 5 ]];then
	z=$((RANDOM%op_sz))
	echo "${opcodes[$z]}"
    else
	z=$((RANDOM%sz2))
	y=$(echo "${expressions[$z]}" |sed 's/[^X]//g' | awk '{print length;}')
	str=$(for c in "${opcodes[@]}";do echo $c;done|shuf -n $y|tr '\n' ','|sed 's/,$//')
	echo $(perl -se 'foreach $x (split(/,/,$b,-1)){$a=~s/X/$x/;}print $a."\n";' -- -a="${expressions[$z]}" -b="$str")
    fi
done|sort|uniq > "${outname1}"

declare -a temp1
declare -a temp2
declare -a temp3
# either ICD9 or ICD10, not both
echo -e "ID\tVISIT\tARRAY\tOPCODE\tICD9\tICD10" > "${outname2}"
for i in $(seq -w 1 $nsamples);do
    n_visits=$((RANDOM%max_visits))
    n_arrays=$((RANDOM%max_arrays))
    p=$((RANDOM%100))
    
    for j in $(seq 0 $n_visits);do
	temp1=()
	temp2=()
	temp3=()
	while read x;do
	    temp1+=($x)
	done < <(echo "$opcodes_str"|tr ',' '\n'|shuf -n $((n_arrays+1)))
	if [[ $p -lt 10 ]];then
	    while read x;do
		temp2+=($x)
	    done < <(echo "$icd9codes_str"|tr ',' '\n'|shuf -n $((n_arrays+1)))
	    while read x;do
		temp3+=($x)
	    done < <(yes "NA"|head -n $((n_arrays+1)))
	else
	    while read x;do
		temp3+=($x)
	    done < <(echo "$icd10codes_str"|tr ',' '\n'|shuf -n $((n_arrays+1)))
	    while read x;do
		temp2+=($x)
	    done < <(yes "NA"|head -n $((n_arrays+1)))
	fi	
	for k in $(seq 0 $n_arrays);do
	    echo $i $j $k "${temp1[$k]}" "${temp2[$k]}" "${temp3[$k]}"
	done
    done
done|tr ' ' '\t' >> "${outname2}"

exit 0
