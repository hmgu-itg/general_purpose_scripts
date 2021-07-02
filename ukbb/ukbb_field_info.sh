#!/usr/bin/env bash

function usage () {
    echo "Get info for a given field"
    echo ""
    exit 0
}

scriptname=$0
scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

declare -A available_projects
declare -A fields
declare -A nacount

config=""
data_path=""
project=""
release=""
infield=""
while getopts "hd:p:r:c:i:" opt; do
    case $opt in
        i)infield=($OPTARG);;
        d)data_path=($OPTARG);;
        c)config=($OPTARG);;
        p)project=($OPTARG);;
        r)release=($OPTARG);;
        h)usage;;
        *)usage;;
    esac
done

if [[ -z "$config" ]];then
    config="${scriptdir}/config.txt"
fi
exitIfEmpty "$release" "ERROR: release (-r) not specified"
exitIfEmpty "$infield" "ERROR: input (-i) not specified"
exitIfNotFile "$config" "ERROR: config file $config does not exist"

readAArray "$config" PROJECTS available_projects
if [[ -z "${available_projects[$project]}" ]];then
    echo "ERROR: project $project is not defined in $config"
    exit 1
fi

if [[ -z "$data_path" ]];then
    readValue "$config" DATA_PATH data_path
fi
exitIfNotDir "$data_path" "ERROR: $data_path is not a directory"
data_path=${data_path%/}

infile="$data_path"/"${available_projects[$project]}"/releases/phenotypes_r"$release".txt.gz
exitIfNotFile "$infile" "ERROR: input file $infile does not exist"
totlines=$(zcat "$infile"|tail -n +3|wc -l)

echo ""
echo "INFO: field: $infield"
echo "INFO: release: $release"
echo "INFO: project: $project"
echo "INFO: config: $config"
echo "INFO: data path: $data_path"
echo "INFO: input file: $infile"
echo ""

getCols "$infile" "zcat" "$infield" fields

if [[ "${#fields[@]}" -eq 0 ]];then
    echo "INFO: $infield not found in $infile header"
    exit 0
fi

echo "INFO: $infield found in columns:" $(join_by , "${!fields[@]}")
echo "INFO: matching fields:" $(join_by , "${fields[@]}")

for c in ${!fields[@]};do
    nacount[$c]=$(zcat "$infile"|cut -f $c|tail -n +3|grep -E "^NA$"|wc -l)
done
echo ""
echo "INFO: NA stats"
for c in ${!fields[@]};do
    echo ${fields[$c]} ${nacount[$c]}/$totlines
done


