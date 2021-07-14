#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
script="${scriptdir}/ukbb_select.sh"

function usage () {
    echo ""
    echo "Script for selecting self-reported OA status from a UKBB release"
    echo ""
    echo "Usage: ukbb_select.sh -r | --release <release>"
    echo "                      -o | --output <output prefix>"
    echo "                      -c | --config <optional: config file; default: config.txt in script directory>"
    echo "                      -h | --help"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

OPTS=$(getopt -o hnc:p:r:o: -l help,names,project:,release:,config:,mean:,majority:,min-missing:,cc:,output: -n 'ukbb_select' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

declare -a ccargs
declare -a majorityargs
declare -a meanargs
declare -a minnaargs

declare -A ccfields
declare -A majorityfields
declare -A meanfields
declare -A minnafields

declare -A available_projects

config=""
release=""
outprefix=""
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
    -r|--release ) release=$2; shift 2 ;;
    -c|--config ) config=$2; shift 2 ;;
    -o|--output ) outprefix=$2; shift 2 ;;
    --) shift ; break ;;
    * ) echo "Unexpected option: $1" ; usage ;;
  esac
done

exitIfEmpty "$release" "ERROR: release not specified"
exitIfEmpty "$outprefix" "ERROR: output prefix not specified"
if [[ -z "$config" ]];then
    config="${scriptdir}"/config.txt
fi

"$script" -p "OA" -r "$release" --cc 20002,1465 -o "$outprefix" -c "$config"
