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
    echo "Preparing input data for FINEMAP"
    echo ""
    echo "Usage: finemap_prepare.sh -i <input.signal.txt> -o <output.prefix> -p <LD.panel.prefix> -a <effect allele column; default:\"allele1\"> -t <threads; default: 1>"
    echo ""
    echo "Input file is space separated, required columns: chromosome"
    echo "                                               : position"
    echo "                                               : rsid"
    echo "                                               : allele1"
    echo "                                               : allele2"
    echo "                                               : n"
    echo "                                               : beta"
    echo "                                               : se"
    echo "                                               : maf"
    exit 0
}

# original ID --> original ID
# new ID --> original ID
# declare -A id_mapping

declare -A input_colnames

# original ID --> effect allele
declare -A ref_alleles

if [[ $# -eq 0 ]];then
    usage
fi

input_fname=""
output_prefix=""
panel_prefix=""
eff_colname="allele1"
threads=1
while getopts "hi:o:p:a:t:" opt; do
    case $opt in
        i)input_fname=($OPTARG);;
        o)output_prefix=($OPTARG);;
        p)panel_prefix=($OPTARG);;
        a)eff_colname=($OPTARG);;
        t)threads=($OPTARG);;
        h)usage;;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

output_prefix=${output_prefix%/}

exitIfEmpty "$input_fname" "ERROR: no input specified"
exitIfEmpty "$output_prefix" "ERROR: no output prefix specified"
exitIfEmpty "$panel_prefix" "ERROR: no LD panel prefix specified"
exitUnlessExists "$input_fname" "ERROR: $input_fname does not exist"
exitUnlessExists "${panel_prefix}".bim "ERROR: ${panel_prefix}.bim does not exist"
exitUnlessExists "${panel_prefix}".bed "ERROR: ${panel_prefix}.bed does not exist"
exitUnlessExists "${panel_prefix}".fam "ERROR: ${panel_prefix}.fam does not exist"

getColNumbersSep "$input_fname" "cat" " " input_colnames

# for i in "${!input_colnames[@]}"
# do
#   echo "name   : $i"
#   echo "column : ${input_colnames[$i]}"
# done

# get column numbers

for n in "chromosome" "position" "rsid" "allele1" "allele2" "beta" "se" "n" "maf";do
    if [[ -v "input_colnames[$n]" ]] ; then
	eval "${n}_column=${input_colnames[$n]}"
    else
	echo "ERROR: could not find $n in the input header"
	exit 1
    fi
done

for n in "chromosome" "position" "rsid" "allele1" "allele2" "beta" "se" "n" "maf";do
    v="${n}_column"
    echo "$n : ${!v}"
done

# exit 0

if [[ -v "input_colnames[$eff_colname]" ]] ; then
    eff_column=${input_colnames[$eff_colname]}
else
    echo "ERROR: could not find $eff_colname in the input header"
    exit 1
fi

# check if input has any duplicate IDs
ndup=$(tail -n +2 $input_fname| cut -d ' ' -f $rsid_column| sort|uniq -d| wc -l)
if [[ $ndup -ne "0" ]];then
    echo "ERROR: there are ID duplicates in input"
    exit 1
fi

# create ID mapping
# while read id;do
#     id_mapping[$id]=$id
#     id_mapping[$(switchIdAlleles $id)]=$id
# done < <(tail -n +2 $input_fname| cut -d ' ' -f $rsid_column)

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
bname=$(basename $(realpath $output_prefix))
tmpdir=$(mktemp -d -p $(dirname $(realpath $output_prefix)) ${bname}_XXXXXXXX)
echo "DEBUG: created $tmpdir"

# effect allele as reference allele
ref_file=${tmpdir}/ref
# for id in "${!id_mapping[@]}"
# do
#     oldID=${id_mapping[${id}]}
#     a=${ref_alleles[$oldID]}
#     echo $id $a
# done > ${ref_file}
for id in "${!ref_alleles[@]}"
do
    a=${ref_alleles[$id]}
    echo $id $a
done > ${ref_file}

# IDs to extract
ex_file=${tmpdir}/ex
# for i in "${!id_mapping[@]}"
# do
#   echo "$i"
# done > ${ex_file}
for i in "${!ref_alleles[@]}"
do
  echo "$i"
done > ${ex_file}

echo "Extracting $(cat ${ex_file}| wc -l) variants"

# calling PLINK
res=${tmpdir}/plink_output
plink --bfile "${panel_prefix}" --extract "${ex_file}" --out "${res}" --reference-allele "${ref_file}"  --r square spaces --write-snplist --threads "$threads"

output_prefix=$(realpath $output_prefix)
ld_fname="${output_prefix}".ld
cp "${res}".ld "$ld_fname"

# create .z file: select lines from the input file that correspond to the snps in plink_output.snplist, in the same order
output_z="${output_prefix}".z
echo  "rsid chromosome position allele1 allele2 maf beta se" > "$output_z"
while read id;do
    awk -v x="$id" -v r="$rsid_column" -v c="$chromosome_column" -v p="$position_column" -v a="$allele1_column" -v d="$allele2_column" -v m="$maf_column" -v b="$beta_column" -v s="$se_column" '{if ($r==x){print x,$c,$p,$a,$d,$m,$b,$s;}}' "$input_fname" >>  "$output_z"
done < <(cat ${res}.snplist)

# N samples for the master file: the max number of samples in the input file
n_samples=$(tail -n +2 "$input_fname"|cut -d ' ' -f $n_column|sort -n|tail -n 1)

# create master file
master_fname="${output_prefix}".master
log_fname="${output_prefix}".log
snp_fname="${output_prefix}".snp
config_fname="${output_prefix}".config
cred_fname="${output_prefix}".cred
echo "z;ld;snp;config;cred;log;n_samples" > "$master_fname"
echo "${output_z};${ld_fname};${snp_fname};${config_fname};${cred_fname};${log_fname};${n_samples}" >> "$master_fname"

echo "DEBUG: deleting $tmpdir"
rm -rf "$tmpdir"

exit 0
