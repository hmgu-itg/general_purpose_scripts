#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
script="${scriptdir}/ukbb_select.sh"
fullsname=$(realpath $0)

function usage () {
    echo ""
    echo "Wrapper script for selecting self reported OA hip/knee/both replacement cases from a UKBB release"
    echo ""
    echo "Usage: OA_select_SR_replacement.sh -r | --release <release of the MAIN dataset>"
    echo "                                   -o | --output <output prefix>"
    echo "                                   -k | --key <one of: TKR/THR/TJR>"
    echo "                                   -c | --config <optional: config file; default: config.txt in script directory>"
    echo "                                   -h | --help"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

OPTS=$(getopt -o hc:r:o:k: -l help,release:,config:,output:,key: -n 'OA_select_SR' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

declare -A keys=( [THR]=1 [TKR]=1 [TJR]=1 )
declare -A keys2=( [THR]=1318 [TKR]=1319 )

config=""
release=""
outprefix=""
key=""
outfile=""
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
    -r|--release ) release=$2; shift 2 ;;
    -k|--key ) key=$2; shift 2 ;;
    -c|--config ) config=$2; shift 2 ;;
    -o|--output ) outprefix=$2; shift 2 ;;
    --) shift ; break ;;
    * ) echo "Unexpected option: $1" ; usage ;;
  esac
done

exitIfEmpty "$key" "ERROR: key not specified"
exitIfEmpty "${keys[$key]}" "ERROR: invalid key $key"
exitIfEmpty "$release" "ERROR: release not specified"
if [[ -z "$config" ]];then
    config="${scriptdir}"/config.txt
fi

# --------------------------------------------------------------------------------------------------------

if [[ ! -z "$outprefix" ]];then
    logfile="$outprefix".log
    outfile="$outprefix".txt.gz
    exitIfExists "$outfile" "ERROR: output file $outfile already exists"
    : > "$logfile"
    sfx="|tee -a $logfile"
else
    sfx="1>&2"
fi

eval echo "Wrapper script: current dir: ${PWD}" "$sfx"
eval echo "Wrapper script: command line: $scriptname ${args[@]}" "$sfx"
eval echo "" "$sfx"

# --------------------------------------------------------------------------------------------------------

if [[ "$key" == "TJR" ]];then
    if [[ -z "${outfile}" ]];then
	cat <("$fullsname" -c "$config" -r "$release" -k "TKR" 2>/dev/null) <("$fullsname" -c "$config" -r "$release" -k "THR" 2>/dev/null)|sort|uniq
    else
	cat <("$fullsname" -c "$config" -r "$release" -k "TKR" 2>/dev/null) <("$fullsname" -c "$config" -r "$release" -k "THR" 2>/dev/null)|sort|uniq|gzip - -c >"$outfile" 
    fi
else
    if [[ -z "${outfile}" ]];then 
	"$script" -p "OA" -r "$release" --cc 20002,1465 --cc 20004,"${keys2[$key]}" -c "$config" | tail -n +2 | awk 'BEGIN{FS="\t";}$2=="1" && $3=="1"{print $1;}'
    else
	"$script" -p "OA" -r "$release" --cc 20002,1465 --cc 20004,"${keys2[$key]}" -c "$config" | tail -n +2 | awk 'BEGIN{FS="\t";}$2=="1" && $3=="1"{print $1;}'| gzip - -c >"$outfile" 
    fi
fi

exit 0
