#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

totalcpus=$(totalCPUs)
echo "Total CPUs: $totalcpus"
#frac=0.5
#echo "Fraction: $frac"
#echo $(getFreeCPUs $frac)
#exit 0

declare -A fn
fn["key1"]="temp_test_XXXXXX"
fn["key2"]="temp_test_XXXXXX"
fn["key3"]="temp_test_XXXXXX"

if createTempFiles "$scriptdir" fn;then
    echo "SUCCESS"
else
    echo "FAILURE"
fi

declare -p fn

for k in "${fn[@]}";do
    if [[ -f "$k" ]];then
	rm -v "$k"
    fi
done

createTempFilesN "$scriptdir" 5 fn "abc_XXX.gz"

echo "in test.sh"
declare -p fn

for k in "${fn[@]}";do
    if [[ -f "$k" ]];then
	rm -v "$k"
    fi
done

exit 0
