#!/usr/bin/env bash

scriptname=$0
args=("$@")

main=$1
icdcodes=$2

scriptdir=$(dirname $(readlink -f $0))
upperdir=$(dirname $scriptdir)
source "${upperdir}/functions.sh"
selectscript="${scriptdir}/check.pl"

tmp_icd9=$(mktemp temp_ICD9_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary file"
    exit 1
fi

tmp_icd10=$(mktemp temp_ICD10_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary file"
    rm "$tmp_icd9"
    exit 1
fi

tmp_data9=$(mktemp temp_data9_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary file"
    rm "$tmp_icd9"
    rm "$tmp_icd10"
    exit 1
fi

tmp_data10=$(mktemp temp_data10_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary file"
    rm "$tmp_icd9"
    rm "$tmp_icd10"
    rm "$tmp_data9"
    exit 1
fi

#----------------------------------------------------------------------------------------------------------------

grep -v "^#" "$icdcodes" | awk 'BEGIN{FS="\t";}$1=="TRAIT" && $2=="ICD9"{print $3;}' > "$tmp_icd9"
grep -v "^#" "$icdcodes" | awk 'BEGIN{FS="\t";}$1=="TRAIT" && $2=="ICD10"{print $3;}' > "$tmp_icd10"
cut -f 1,4 "$main"|datamash -s -g 1 collapse 2 > "$tmp_data9"
cut -f 1,5 "$main"|datamash -s -g 1 collapse 2 > "$tmp_data10"

cat <("$selectscript" "$tmp_icd9" "$tmp_data9" 1) <("$selectscript" "$tmp_icd10" "$tmp_data10" 1)|sort|uniq 

rm "$tmp_icd9" "$tmp_icd10" "$tmp_data9" "$tmp_data10"
exit 0
