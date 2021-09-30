#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

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

exit 0
