#!/usr/bin/env bash

function usage () {
    echo "List available UKBB projects and releases"
    echo ""
    echo "Usage: ukbb_release_info.sh -c <config.txt: optional, default: config.txt in this script directory>"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

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
    encoding="${scriptdir}"/config.txt
    echo "INFO: no config file specified, trying to use $config"
fi
exitIfNotFile "$config" "ERROR: config $config does not exist"
readAArray "$config" PROJECTS projects
readValue "$config" DATA_PATH data_path
