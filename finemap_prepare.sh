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
    echo "Usage: finemap_prepare.sh -i <input.signal.txt> -o <output.prefix> -p <LD.panel.prefix> -a <effect allele column; default:\"allele1\">"
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
declare -A id_mapping

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
while getopts "hi:o:p:a:" opt; do
    case $opt in
        i)input_fname=($OPTARG);;
        o)output_prefix=($OPTARG);;
        p)panel_prefix=($OPTARG);;
        a)eff_colname=($OPTARG);;
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

# sep=" "

getColNumbersSep "$input_fname" "cat" " " input_colnames

# for i in "${!input_colnames[@]}"
# do
#   echo "name   : $i"
#   echo "column : ${input_colnames[$i]}"
# done

# get column numbers

if [[ -v "input_colnames[$eff_colname]" ]] ; then
    eff_column=${input_colnames[$eff_colname]}
else
    echo "ERROR: could not find $eff_colname in the input header"
    exit 1
fi

echo "$eff_colname --> $eff_column"

if [[ -v "input_colnames[chromosome]" ]] ; then
    chr_column=${input_colnames["chromosome"]}
else
    echo "ERROR: could not find \"chromosome\" in the input header"
    exit 1
fi

if [[ -v "input_colnames[position]" ]] ; then
    pos_column=${input_colnames["position"]}
else
    echo "ERROR: could not find \"position\" in the input header"
    exit 1
fi

if [[ -v "input_colnames[rsid]" ]] ; then
    id_column=${input_colnames["rsid"]}
else
    echo "ERROR: could not find \"rsid\" in the input header"
    exit 1
fi

if [[ -v "input_colnames[allele1]" ]] ; then
    a1_column=${input_colnames["allele1"]}
else
    echo "ERROR: could not find \"allele1\" in the input header"
    exit 1
fi

if [[ -v "input_colnames[allele2]" ]] ; then
    a2_column=${input_colnames["allele2"]}
else
    echo "ERROR: could not find \"allele2\" in the input header"
    exit 1
fi

if [[ -v "input_colnames[maf]" ]] ; then
    maf_column=${input_colnames["maf"]}
else
    echo "ERROR: could not find \"maf\" in the input header"
    exit 1
fi

if [[ -v "input_colnames[beta]" ]] ; then
    beta_column=${input_colnames["beta"]}
else
    echo "ERROR: could not find \"beta\" in the input header"
    exit 1
fi

if [[ -v "input_colnames[se]" ]] ; then
    se_column=${input_colnames["se"]}
else
    echo "ERROR: could not find \"se\" in the input header"
    exit 1
fi

if [[ -v "input_colnames[n]" ]] ; then
    n_column=${input_colnames["n"]}
else
    echo "ERROR: could not find \"n\" in the input header"
    exit 1
fi

# check if input has any duplicate positions
ndup=$(tail -n +2 $input_fname| cut -d ' ' -f $chr_column,$pos_column| tr ' ' ':'| sort|uniq -d| wc -l)

if [[ $ndup -ne "0" ]];then
    echo "ERROR: there are chr:pos duplicates in input"
    exit 1
fi

# create ID mapping
while read id;do
    id_mapping[$id]=$id
    id_mapping[$(switchIdAlleles $id)]=$id
done < <(tail -n +2 $input_fname| cut -d ' ' -f $id_column)

# create "old ID" --> "ref allele" mapping
if [[ "$id_column" -lt "$eff_column" ]];then
    while read id a;do
	ref_alleles[$id]=$a
    done < <(tail -n +2 $input_fname| cut -d ' ' -f $id_column,$eff_column)
else
    while read a id;do
	ref_alleles[$id]=$a
    done < <(tail -n +2 $input_fname| cut -d ' ' -f $id_column,$eff_column)
fi

id="4:1234567_A_G"
newID=$(switchIdAlleles $id)
echo $id
echo $newID

# create temp dir
tmpdir=$(mktemp -d -p $(dirname $(realpath $output_prefix)) ${output_prefix}_XXXXXXXX)
echo "DEBUG: created $tmpdir"

# for i in "${!ref_alleles[@]}"
# do
#   echo "$i" "${ref_alleles[$i]}"
# done

# effect allele as reference allele
ref_file=${tmpdir}/ref
for id in "${!id_mapping[@]}"
do
    oldID=${id_mapping[${id}]}
    # echo $id $oldID
    a=${ref_alleles[$oldID]}
    echo $id $a
done > ${ref_file}

# IDs to extract
ex_file=${tmpdir}/ex
for i in "${!id_mapping[@]}"
do
  echo "$i"
done > ${ex_file}

# calling PLINK
res=${tmpdir}/plink_output
plink --bfile "${panel_prefix}" --extract "${ex_file}" --out "${res}" --reference-allele "${ref_file}"  --r square spaces --write-snplist

output_prefix=$(realpath $output_prefix)
ld_fname="${output_prefix}".ld
cp "${res}".ld "$ld_fname"

# create .z file: select lines from the input file that correspond to the snps in plink_output.snplist, in the same order
output_z="${output_prefix}".z
echo  "rsid chromosome position allele1 allele2 maf beta se" > "$output_z"
while read id;do
    oldID=${id_mapping[$id]}
    #echo "id=$oldID rs=$id_column chr=$chr_column pos=$pos_column a1=$a1_column a2=$a2_column maf=$maf_column beta=$beta_column se=$se_column"
    awk -v x="$oldID" -v r="$id_column" -v c="$chr_column" -v p="$pos_column" -v a="$a1_column" -v c="$a2_column" -v m="$maf_column" -v b="$beta_column" -v s="$se_column" '{if ($r==x){print x,$c,$p,$a,$c,$m,$b,$s;}}' "$input_fname" >>  "$output_z"
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
