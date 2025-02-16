#!/usr/bin/env bash

function usage () {
    echo "List available UKBB projects and releases"
    echo ""
    echo "Usage: ukbb_release_info.sh -c <config.txt: optional, default: config.txt in this script's directory>"
    echo ""
    exit 0
}

scriptname=$0
scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

declare -A projects

data_path=""
config=""
while getopts "hc:" opt; do
    case $opt in
        c)config=($OPTARG);;
        h)usage;;
        *)usage;;
    esac
done

if [[ -z "$config" ]];then
    config="${scriptdir}"/config.txt
    echo "INFO: no config file specified, using $config"
fi
exitIfNotFile "$config" "ERROR: config $config does not exist"
readAArray "$config" PROJECTS projects
readValue "$config" DATA_PATH data_path
echo -e "PROJECT\tDATASET\tRELEASE\tCREATED\tSAMPLES\tFIELDS"
for p in "${!projects[@]}";do
    dir="$data_path"/"${projects[$p]}"/"releases"
    for f in $(find "$dir" -maxdepth 1 -type f -name "phenotypes_r*.txt.gz" ! -name "*.intersection.txt.gz" | sort);do
	rel=$(getRelease "$f" "zcat")
	cr=$(getCreated "$f" "zcat")
	samples=$(zcat "$f" | tail -n +2 | wc -l)
	fields=$(zcat "$f" | head -n 1 | tr '\t' '\n' | wc -l)
	fields=$((fields-3)) # exclude f.eid, RELEASE, CREATED
	echo -e "${p}\tmain\t${rel}\t${cr}\t${samples}\t${fields}"
    done
    for f in $(find "$dir" -maxdepth 1 -type f -name "hesin_r*.tar.gz"| sort);do
	rel=$(tar -zxf "$f" RELEASE -O)
	cr=$(tar -zxf "$f" CREATED -O)
	echo -e "${p}\thesin\t${rel}\t${cr}\tNA\tNA"
    done
done

exit 0
