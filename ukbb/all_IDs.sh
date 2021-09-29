#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

function usage () {
    echo ""
    echo "Selecting all sample IDs from a UKBB release"
    echo ""
    echo "Usage: all_IDs.sh -p | --project <project name>"
    echo "                  -r | --release <release of the MAIN dataset>"
    echo "                  -o | --output <optional: output prefix; default: STDOUT>"
    echo "                  -c | --config <optional: config file; default: config.txt in script directory>"
    echo "                  -h | --help"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

OPTS=$(getopt -o hc:r:o:p: -l help,project:,release:,config:,output: -n 'all_IDs' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

declare -A available_projects
config=""
project=""
release=""
outprefix=""
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
    -r|--release ) release=$2; shift 2 ;;
    -p|--project ) project=$2; shift 2 ;;
    -c|--config ) config=$2; shift 2 ;;
    -o|--output ) outprefix=$2; shift 2 ;;
    --) shift ; break ;;
    * ) echo "Unexpected option: $1" ; usage ;;
  esac
done

exitIfEmpty "$project" "ERROR: project name not specified"
exitIfEmpty "$release" "ERROR: MAIN release not specified"
if [[ -z "$config" ]];then
    config="${scriptdir}"/config.txt
fi
readValue "$config" DATA_PATH data_path
exitIfNotFile "$config" "ERROR: config $config does not exist"
readAArray "$config" PROJECTS available_projects
if [[ -z "${available_projects[$project]}" ]];then
    echo "ERROR: project $project is not defined in $config"
    exit 1
fi
infile="$data_path"/"${available_projects[$project]}"/releases/phenotypes_r"$release".txt.gz
exitIfNotFile "$infile" "ERROR: project file $infile does not exist"

if [[ ! -z "$outprefix" ]];then
    logfile="$outprefix".log
    outfile="$outprefix".txt.gz
    exitIfExists "$outfile" "ERROR: output file $outfile already exists"
    : > "$logfile"
    sfx="|tee -a $logfile"
else
    sfx="1>&2"
fi

eval date "+%d-%b-%Y:%H-%M-%S" "$sfx"
eval echo "" "$sfx"
eval echo "Current dir: ${PWD}" "$sfx"
eval echo "Command line: $scriptname ${args[@]}" "$sfx"
eval echo "" "$sfx"

eval echo "INFO: project: $project" "$sfx"
eval echo "INFO: release: $release" "$sfx"
eval echo "INFO: project file: $infile" "$sfx"
eval echo "INFO: output prefix: $outprefix" "$sfx"
eval echo "INFO: output file: $outfile" "$sfx"
eval echo "" "$sfx"

eid_col=$(getColNum "$infile" "f.eid" "zcat")
if [[ -z "$outfile" ]];then
    zcat "$infile"| tail -n +3|cut -f "$eid_col"
else
    zcat "$infile"| tail -n +3|cut -f "$eid_col" | gzip - -c > "$outfile"
fi

exit 0
