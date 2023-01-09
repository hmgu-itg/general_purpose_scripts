#!/usr/bin/bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $scriptname))
source "${scriptdir}/ukbb/functions.sh"

# id: ref. panel ID: "chr:pos:A1:A2"
# ea: EA from the signal input file
function getEffectAllele {
    local id=$1
    local ea=$2
    local nea=$3
    
    echo $id|tr '_' ':'|perl -slne '@a=split(/:/);$a1=$a[2];$a2=$a[3];if ($x eq $a1 || $x eq $a2){print $x;exit 0;} if ($x eq "D"){if (length($a1)>length($a2)){print $a2;}elsif(length($a2)>length($a1)){print $a1;}else{print "NA";}exit 0;} if ($x eq "I"){if (length($a1)>length($a2)){print $a1;}elsif(length($a2)>length($a1)){print $a2;}else{print "NA";}exit 0;} if (length($x)>length($z)){if ($a1 eq "I" || $a2 eq "I"){print "I";}else{print "NA";}exit 0;} if (length($x)<length($z)){if ($a1 eq "D" || $a2 eq "D"){print "D";}else{print "NA";}exit 0;} print "NA";' -- -x=$ea -z=$nea
}

function usage () {
    echo ""
    echo "Preparing input data for SUSIE"
    echo ""
    echo "Usage: susie_prepare.sh -i <input.signal.txt> -o <output.prefix> -p <LD.panel.prefix> -m <ID.mapping.txt; optional> -a <effect allele column; default:\"allele1\"> -n <non-effect allele column; default:\"allele2\"> -t <threads; default: 1> -k <keep temp files> -c <MAC; default: 0>"
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

#ref. panel ID --> input ID
declare -A id_mapping
declare -A input_colnames
# input ID --> effect allele
declare -A ref_alleles
# input ID --> non-effect allele
declare -A nref_alleles

if [[ $# -eq 0 ]];then
    usage
fi

mapping_fname=""
input_fname=""
output_prefix=""
panel_prefix=""
eff_colname="allele1"
neff_colname="allele2"
keep="NO"
threads=1
mac=0
while getopts "hi:t:o:p:a:m:kc:" opt; do
    case $opt in
        i)input_fname=($OPTARG);;
        m)mapping_fname=($OPTARG);;
        o)output_prefix=($OPTARG);;
        p)panel_prefix=($OPTARG);;
        a)eff_colname=($OPTARG);;
        n)neff_colname=($OPTARG);;
        t)threads=($OPTARG);;
        c)mac=($OPTARG);;
	k)keep="YES";;
        h)usage;;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

exitIfEmpty "$input_fname" "ERROR: no input specified"
# exitIfEmpty "$mapping_fname" "ERROR: no ID mapping specified"
exitIfEmpty "$output_prefix" "ERROR: no output prefix specified"
exitIfEmpty "$panel_prefix" "ERROR: no LD panel prefix specified"
exitUnlessExists "$input_fname" "ERROR: $input_fname does not exist"
# exitUnlessExists "$mapping_fname" "ERROR: $mapping_fname does not exist"
exitUnlessExists "${panel_prefix}".bim "ERROR: ${panel_prefix}.bim does not exist"
exitUnlessExists "${panel_prefix}".bed "ERROR: ${panel_prefix}.bed does not exist"
exitUnlessExists "${panel_prefix}".fam "ERROR: ${panel_prefix}.fam does not exist"

getColNumbersSep "$input_fname" "cat" " " input_colnames

log_fname="${output_prefix}".log

: > "${log_fname}"

echo "INPUT: $input_fname" | tee -a "${log_fname}"
echo "REF PANEL PREFIX: $panel_prefix" | tee -a "${log_fname}"
echo "ID MAPPING: $mapping_fname" | tee -a "${log_fname}"
echo "OUTPUT: $output_prefix" | tee -a "${log_fname}"
echo "KEEP TEMP FILES: $keep" | tee -a "${log_fname}"
echo "MAC THRESHOLD: $mac" | tee -a "${log_fname}"
echo "THREADS: $threads" | tee -a "${log_fname}"
echo "EA COLNAME: $eff_colname" | tee -a "${log_fname}"
echo "NEA COLNAME: $neff_colname" | tee -a "${log_fname}"

# get column numbers
for n in "chromosome" "position" "rsid" "allele1" "allele2" "beta" "se" "n";do
    if [[ -v "input_colnames[$n]" ]] ; then
	eval "${n}_column=${input_colnames[$n]}"
    else
	echo "ERROR: could not find $n in the input header" | tee -a "${log_fname}"
	exit 1
    fi
done

if [[ -v "input_colnames[$eff_colname]" ]] ; then
    eff_column=${input_colnames[$eff_colname]}
else
    echo "ERROR: could not find $eff_colname in the input header" | tee -a "${log_fname}"
    exit 1
fi

if [[ -v "input_colnames[$neff_colname]" ]] ; then
    neff_column=${input_colnames[$neff_colname]}
else
    echo "ERROR: could not find $neff_colname in the input header" | tee -a "${log_fname}"
    exit 1
fi

echo "EA COLUMN: $eff_column" | tee -a "${log_fname}"
echo "NEA COLUMN: $neff_column" | tee -a "${log_fname}"

# total number of input variants
echo "Input variants: $(tail -n +2 ${input_fname}|wc -l)" | tee -a "${log_fname}"

# ID mapping
# ref. panel ID --> input ID
# identity mapping for variants common to input and ref. panel, if no mapping_fname is provided
if [[ -z "${mapping_fname}" ]];then
    while read id;do
	id_mapping[$id]=$id
    done < <(cat <(tail -n +2 $input_fname|cut -d ' ' -f $rsid_column) <(sed 's/\t/ /g' "${panel_prefix}".bim|cut -d ' ' -f 2)|sort|uniq -d)
else
    while read id id2;do
	id_mapping[$id2]=$id
    done < <(join -1 1 -2 1 <(tail -n +2 $input_fname|cut -d ' ' -f $rsid_column|sort) <(grep -v "^#" $mapping_fname|sort -t ' ' -k1,1))
fi

echo "Mapped IDs in input: ${#id_mapping[@]}" | tee -a "${log_fname}"
echo $(date) "DEBUG: ID mapping done" | tee -a "${log_fname}"

# "input ID" --> "effect allele" mapping
if [[ "$rsid_column" -lt "$eff_column" ]];then
    while read id a;do
	ref_alleles[$id]=$a
    done < <(tail -n +2 $input_fname| cut -d ' ' -f $rsid_column,$eff_column)
else
    while read a id;do
	ref_alleles[$id]=$a
    done < <(tail -n +2 $input_fname| cut -d ' ' -f $rsid_column,$eff_column)
fi
echo $(date) "DEBUG: EA mapping done" | tee -a "${log_fname}"

# "input ID" --> "non-effect allele" mapping
if [[ "$rsid_column" -lt "$neff_column" ]];then
    while read id a;do
	nref_alleles[$id]=$a
    done < <(tail -n +2 $input_fname| cut -d ' ' -f $rsid_column,$neff_column)
else
    while read a id;do
	nref_alleles[$id]=$a
    done < <(tail -n +2 $input_fname| cut -d ' ' -f $rsid_column,$neff_column)
fi
echo $(date) "DEBUG: NEA mapping done" | tee -a "${log_fname}"

# create temp dir
output_fname="${output_prefix}".txt
bname=$(basename $(realpath $output_fname))
tmpdir=$(mktemp -d -p $(dirname $(realpath $output_fname)) ${bname}_XXXXXXXX)

# use effect allele as PLINK reference allele
ref_file=${tmpdir}/ref
# IDs to extract
ex_file=${tmpdir}/ex
for id in "${!id_mapping[@]}"
do
    id2=${id_mapping[${id}]}
    a=${ref_alleles[${id2}]}
    b=${nref_alleles[${id2}]}
    a2=$(getEffectAllele $id $a $b)
    if [[ $a2 == "NA" ]];then
	echo "ERROR: could not get effect allele for $id ($a); skipping" | tee -a "${log_fname}"
    else
	echo $id $a2 >> ${ref_file}
	echo $id >> ${ex_file}
    fi
done

echo $(date) "DEBUG: reference alleles saved in ${ref_file}" | tee -a "${log_fname}"
echo $(date) "DEBUG: IDs to be extracted saved in ${ex_file}" | tee -a "${log_fname}"

# using PLINK to create correlation matrix
res=${tmpdir}/plink_output
plink --threads $threads --bfile "${panel_prefix}" --extract "${ex_file}" --out "${res}" --reference-allele "${ref_file}"  --r square spaces --write-snplist --mac $mac 1>${res}.log 2>&1

if [[ $? -ne 0 ]];then
    echo "ERROR: PLINK failed" | tee -a "${log_fname}"
    exit 1
else
    echo $(date) "DEBUG: PLINK finished" | tee -a "${log_fname}"
fi

# number of variants in PLINK output
n_out=$(cat ${res}.snplist| wc -l)
echo "${n_out} variants in PLINK output" | tee -a "${log_fname}"

# temp file for the correlation matrix, with header
c_file=${tmpdir}/corr
perl -sle 'print join(" ",map { "C".$_ } (1 .. $n));' -- -n=${n_out} > "${c_file}"
cat "${res}".ld >> "${c_file}"

# N samples: the same for each SNP, max number of samples in the input file
n_samples=$(tail -n +2 "$input_fname"|cut -d ' ' -f $n_column|sort -n|tail -n 1)

# select lines from the input file that correspond to the snps in plink_output.snplist, in the same order, save in the temp file
t_file=${tmpdir}/part1
echo  "ID N beta se" > "${t_file}" | tee -a "${log_fname}"
while read id;do
    oldID=${id_mapping[$id]}
    awk -v x="$oldID" -v r="$rsid_column" -v n="$n_samples" -v b="$beta_column" -v s="$se_column" '{if ($r==x){print x,n,$b,$s;exit 0;}}' "$input_fname" >>  "${t_file}"
done < <(cat ${res}.snplist)
echo $(date) "DEBUG: output part 1 finished" | tee -a "${log_fname}"

n2=$(tail -n +2 ${t_file}|wc -l)
if [[ "${n_out}" -ne "${n2}" ]];then
    echo "ERROR: could not match all ${n_out} PLINK variants (matched variants: ${n2})" | tee -a "${log_fname}"
    exit 1
fi

# report unmatched variant entries
while read id;do
    awk -v x="$id" -v r="$rsid_column" 'BEGIN{OFS="\t";}{if ($r==x){print "NOT INCLUDED",$0;exit 0;}}' "$input_fname" >> "${log_fname}"
done < <(cat <(tail -n +2 ${t_file}|cut -d ' ' -f 1) <(tail -n +2 ${input_fname}|cut -d ' ' -f 1)|sort|uniq -u)

# combine two temp files
paste -d ' ' "${t_file}" "${c_file}" > "${output_fname}"
echo $(date) "DEBUG: output finished; compressing" | tee -a "${log_fname}"
pigz -f -p $threads "${output_fname}"
echo $(date) "DEBUG: compressing done" | tee -a "${log_fname}"

#echo "DEBUG: deleting ${tmpdir}"
if [[ $keep == "NO" ]];then
    rm -rf "${tmpdir}"
fi

exit 0
