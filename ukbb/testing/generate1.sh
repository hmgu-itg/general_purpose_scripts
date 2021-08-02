#!/usr/bin/env bash

scriptdir=$(dirname $(readlink -f $0))
upperdir=$(dirname $scriptdir)
source "${upperdir}/functions.sh"
scriptname="${scriptdir}/process_chunk.pl"

nsamples=$1
outname1=$2
outname2=$3
outname3=$4

declare -ir max_arrays=10
declare -ir max_visits=10
declare -r p1=0.5
declare -a opcodes
declare -a icd10codes
declare -a icd9codes
declare -ir out_opcodes=50
declare -ir out_icd9codes=20
declare -ir out_icd10codes=100
declare -ir ncases=200
declare -ir n2cases=50

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
icd9_sz="${#icd9codes[@]}"
icd10_sz="${#icd10codes[@]}"
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
done < <(echo "$icd9codes_str"|tr ',' '\n'|shuf -n $((out_icd10codes)))
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
done > "${outname1}"
echo "expressions done"
date "+%d-%b-%Y:%H-%M-%S"

# echo "CASES:" $(join_by , "${!selected_cases[@]}")
# echo "CASES2:" $(join_by , "${!selected2_cases[@]}")
# echo "ICD9:" $(join_by , "${!selected_icd9[@]}")
# echo "ICD10:" $(join_by , "${!selected_icd10[@]}")
# echo "OPCODES:" $(join_by , "${!selected_opcodes[@]}")

opc=$(echo $(join_by , "${!selected_opcodes[@]}")|tr ' ()' '_:?')

seq -w 1 $nsamples | parallel --pipe --block 1M -N1000 "$scriptname" $(join_by , "${!selected_cases[@]}") $(join_by , "${!selected2_cases[@]}") $(join_by , "${!selected_icd9[@]}") $(join_by , "${!selected_icd10[@]}") "$opc" > "$outname2"

echo "main done"
date "+%d-%b-%Y:%H-%M-%S"

exit 0
