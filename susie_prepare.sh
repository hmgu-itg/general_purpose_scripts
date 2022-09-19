#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $scriptname))
source "${scriptdir}/ukbb/functions.sh"

function switchIdAlleles {
    local id=$1
    
    echo $id|perl -lne '/^(.*)_([[:alpha:]])_([[:alpha:]])$/;print $1."_".$3."_".$2;'
}

function usage () {
    echo ""
    echo "Preparing input data for SUSIE"
    echo ""
    echo "Usage: susie_prepare.sh -i <input.signal.txt> -o <output> -p <LD.panel.prefix> -a <effect allele column; default:\"allele1\">"
    echo ""
    echo "Input file is space separated, required columns: chromosome"
    echo "                                               : position"
    echo "                                               : rsid"
    echo "                                               : allele1"
    echo "                                               : allele2"
    echo "                                               : n"
    echo "                                               : beta"
    echo "                                               : se"
    exit 0
}

# original ID --> original ID
# new ID --> original ID
declare -A id_mapping

declare -A input_colnames

# original ID --> effect allele
declare -A ref_alleles

if [[ $# -eq 0 ]];then
    usage
fi

input_fname=""
output_fname=""
panel_prefix=""
eff_colname="allele1"
while getopts "hi:o:p:a:" opt; do
    case $opt in
        i)input_fname=($OPTARG);;
        o)output_fname=($OPTARG);;
        p)panel_prefix=($OPTARG);;
        a)eff_colname=($OPTARG);;
        h)usage;;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

exitIfEmpty "$input_fname" "ERROR: no input specified"
exitIfEmpty "$output_fname" "ERROR: no output prefix specified"
exitIfEmpty "$panel_prefix" "ERROR: no LD panel prefix specified"
exitUnlessExists "$input_fname" "ERROR: $input_fname does not exist"
exitUnlessExists "${panel_prefix}".bim "ERROR: ${panel_prefix}.bim does not exist"
exitUnlessExists "${panel_prefix}".bed "ERROR: ${panel_prefix}.bed does not exist"
exitUnlessExists "${panel_prefix}".fam "ERROR: ${panel_prefix}.fam does not exist"

getColNumbersSep "$input_fname" "cat" " " input_colnames

echo "INPUT: $input_fname"

# get column numbers
for n in "chromosome" "position" "rsid" "allele1" "allele2" "beta" "se" "n";do
    if [[ -v "input_colnames[$n]" ]] ; then
	eval "${n}_column=${input_colnames[$n]}"
    else
	echo "ERROR: could not find $n in the input header"
	exit 1
    fi
done

if [[ -v "input_colnames[$eff_colname]" ]] ; then
    eff_column=${input_colnames[$eff_colname]}
else
    echo "ERROR: could not find $eff_colname in the input header"
    exit 1
fi

# total number of input variants
echo "Input variants: $(tail -n +2 ${input_fname}|wc -l)"

# check if input has any duplicate positions
ndup=$(tail -n +2 $input_fname| cut -d ' ' -f $chromosome_column,$position_column| tr ' ' ':'| sort|uniq -d| wc -l)

if [[ $ndup -ne "0" ]];then
    echo "ERROR: there are chr:pos duplicates in input"
    exit 1
fi

# create ID mapping
while read id;do
    id_mapping[$id]=$id
    id_mapping[$(switchIdAlleles $id)]=$id
done < <(tail -n +2 $input_fname| cut -d ' ' -f $rsid_column)

# create "old ID" --> "ref allele" mapping
if [[ "$rsid_column" -lt "$eff_column" ]];then
    while read id a;do
	ref_alleles[$id]=$a
    done < <(tail -n +2 $input_fname| cut -d ' ' -f $rsid_column,$eff_column)
else
    while read a id;do
	ref_alleles[$id]=$a
    done < <(tail -n +2 $input_fname| cut -d ' ' -f $rsid_column,$eff_column)
fi

# create temp dir
bname=$(basename $(realpath $output_fname))
tmpdir=$(mktemp -d -p $(dirname $(realpath $output_fname)) ${bname}_XXXXXXXX)
#echo "DEBUG: created $tmpdir"

# use effect allele as reference allele
ref_file=${tmpdir}/ref
for id in "${!id_mapping[@]}"
do
    oldID=${id_mapping[${id}]}
    a=${ref_alleles[$oldID]}
    echo $id $a
done > ${ref_file}

# IDs to extract
ex_file=${tmpdir}/ex
for i in "${!id_mapping[@]}"
do
  echo "$i"
done > ${ex_file}

#echo "Extracting $(cat ${ex_file}| wc -l) variants (including IDs with switched alleles)"

# using PLINK to create correlation matrix
res=${tmpdir}/plink_output
plink --bfile "${panel_prefix}" --extract "${ex_file}" --out "${res}" --reference-allele "${ref_file}"  --r square spaces --write-snplist 1>/dev/null 2>&1

if [[ $? -ne 0 ]];then
    echo "ERROR: PLINK failed"
    exit 1
fi

# number of variants in PLINK output
n_out=$(cat ${res}.snplist| wc -l)
echo "${n_out} variants in PLINK output"

# temp file for the correlation matrix, with header
c_file=${tmpdir}/corr
perl -sle 'print join(" ",map { "C".$_ } (1 .. $n));' -- -n=${n_out} > "${c_file}"
cat "${res}".ld >> "${c_file}"

# N samples: the same for each SNP, max number of samples in the input file
n_samples=$(tail -n +2 "$input_fname"|cut -d ' ' -f $n_column|sort -n|tail -n 1)

# select lines from the input file that correspond to the snps in plink_output.snplist, in the same order, save in the temp file
t_file=${tmpdir}/part1
echo  "ID N beta se" > "${t_file}"
while read id;do
    oldID=${id_mapping[$id]}
    awk -v x="$oldID" -v r="$rsid_column" -v n="$n_samples" -v b="$beta_column" -v s="$se_column" '{if ($r==x){print x,n,$b,$s;exit 0;}}' "$input_fname" >>  "${t_file}"
done < <(cat ${res}.snplist)

n2=$(tail -n +2 ${t_file}|wc -l)
if [[ "${n_out}" -ne "${n2}" ]];then
    echo "ERROR: could not match all ${n_out} PLINK variants (matched variants: ${n2})"
    exit 1
fi

# report unmatched variant entries
while read id;do
    awk -v x="$id" -v r="$rsid_column" 'BEGIN{OFS="\t";}{if ($r==x){print "NOT INCLUDED",$0;exit 0;}}' "$input_fname"
done < <(cat <(tail -n +2 ${t_file}|cut -d ' ' -f 1) <(tail -n +2 ${input_fname}|cut -d ' ' -f 1)|sort|uniq -u)

# combine two temp files
paste -d ' ' "${t_file}" "${c_file}" > "${output_fname}"
pigz -f -p 24 "${output_fname}"

#echo "DEBUG: deleting ${tmpdir}"
rm -rf "${tmpdir}"

exit 0
