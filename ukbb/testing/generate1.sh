#!/usr/bin/env bash

scriptdir=$(dirname $(readlink -f $0))
upperdir=$(dirname $scriptdir)
source "${upperdir}/functions.sh"
scriptname="${scriptdir}/process_chunk.pl"

nsamples=$1
ncases=$2
n2cases=$3
outname1=$4
outname2=$5
outname3=$6

declare -a opcodes
declare -a icd10codes
declare -a icd9codes
declare -ir out_opcodes=200
declare -ir out_icd9codes=50
declare -ir out_icd10codes=200

declare -a expressions=("X And (X Or X Or X Or (X AND X) Or X)" "X And (X Or X Or X)")

# ------------------------------------------------------------------------------------------------------------------

for c in {A..B};do
    for i in $(seq -w 0 999);do
	opcodes+=($c$i)
    done
done

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

declare -A selected_cases
declare -A selected2_cases
while read s;do
    selected_cases[$s]=1
done < <(seq -w 1 $nsamples|shuf -n $ncases)
while read s;do
    selected2_cases[$s]=1
done < <(echo $(join_by , "${!selected_cases[@]}")|tr ',' '\n'|shuf -n $n2cases)

date "+%d-%b-%Y:%H-%M-%S"
# ICD codes
declare -A selected_icd9
declare -A selected_icd10
: > "$outname3"
while read x;do
    selected_icd10["$x"]=1
    echo "TRAIT" "ICD10" "$x"| tr ' ' '\t' >> "$outname3"
done < <(echo "$icd10codes_str"|tr ',' '\n'|shuf -n $((out_icd10codes)))
while read x;do
    selected_icd9["$x"]=1
    echo "TRAIT" "ICD9" "$x"| tr ' ' '\t' >> "$outname3"
done < <(echo "$icd9codes_str"|tr ',' '\n'|shuf -n $((out_icd9codes)))
echo "ICD done"
date "+%d-%b-%Y:%H-%M-%S"

declare -A selected_opcodes
# opcode expressions
for (( i=0; i<$out_opcodes; i++ ));do
    x=$((RANDOM%10))
    if [[ $x -lt 5 ]];then
	z=$((RANDOM%op_sz))
	if [[ -z "${selected_opcodes[${opcodes[$z]}]}" ]];then
	    selected_opcodes["${opcodes[$z]}"]=1
	fi
    else
	z=$((RANDOM%sz2))
	y=$(echo "${expressions[$z]}" |sed 's/[^X]//g' | awk '{print length;}')
	str=$(for c in "${opcodes[@]}";do echo $c;done|shuf -n $y|tr '\n' ','|sed 's/,$//')
	s=$(perl -se 'foreach $x (split(/,/,$b,-1)){$a=~s/X/$x/;}print $a."\n";' -- -a="${expressions[$z]}" -b="$str")
	if [[ -z "${selected_opcodes[$s]}" ]];then
	    selected_opcodes["$s"]=1
	fi
    fi
done
for x in "${!selected_opcodes[@]}";do
    echo $x
done|sort|uniq > "${outname1}"
echo "expressions done"
date "+%d-%b-%Y:%H-%M-%S"

opc=$(echo $(join_by , "${!selected_opcodes[@]}")|tr ' ()' '_:?')

seq -w 1 $nsamples | parallel --pipe --block 1M -N1000 "$scriptname" $(join_by , "${!selected_cases[@]}") $(join_by , "${!selected2_cases[@]}") $(join_by , "${!selected_icd9[@]}") $(join_by , "${!selected_icd10[@]}") "$opc" > "$outname2"

echo "main done"
date "+%d-%b-%Y:%H-%M-%S"

exit 0
