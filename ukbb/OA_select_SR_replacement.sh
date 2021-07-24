#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
script="${scriptdir}/ukbb_select.sh"

function usage () {
    echo ""
    echo "Wrapper script for selecting self reported OA hip/knee/both replacement cases from a UKBB release"
    echo ""
    echo "Usage: ukbb_select.sh -r | --release <release of the MAIN dataset>"
    echo "                      -o | --output <output prefix>"
    echo "                      -k | --key <one of: TKR/THR/TJR>"
    echo "                      -c | --config <optional: config file; default: config.txt in script directory>"
    echo "                      -h | --help"
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

config=""
release=""
outprefix=""
key=""
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
exitIfEmpty "${keys[$key]}" "ERROR: wrong key $key"
exitIfEmpty "$release" "ERROR: release not specified"
exitIfEmpty "$outprefix" "ERROR: output prefix not specified"
if [[ -z "$config" ]];then
    config="${scriptdir}"/config.txt
fi

outfile="${outprefix}".txt.gz
exitIfExists "$outfile" "ERROR: output file $outfile already exists"
logfile="${outprefix}".log

: > "$logfile"
echo "Wrapper script: current dir: ${PWD}"|tee -a "$logfile"
echo "Wrapper script: command line: $scriptname ${args[@]}"|tee -a "$logfile"
echo ""|tee -a "$logfile"

case $key in
    THR)
	"$script" -p "OA" -r "$release" --cc 20002,1465 --cc 20004,1318 -c "$config" 2>>"$logfile" | tail -n +2 | awk 'BEGIN{FS="\t";}$2=="1" && $3=="1"{print $1;}'| gzip - -c >"$outfile" 
	;;
    TKR)
	"$script" -p "OA" -r "$release" --cc 20002,1465 --cc 20004,1319 -c "$config" 2>>"$logfile" | tail -n +2 | awk 'BEGIN{FS="\t";}$2=="1" && $3=="1"{print $1;}'| gzip - -c >"$outfile" 
	;;
    TJR)
	"$script" -p "OA" -r "$release" --cc 20002,1465 --cc 20004,1318 --cc 20004,1319 -c "$config" 2>>"$logfile" | awk 'BEGIN{FS="\t";}{if (NR==1){for (i=1;i<=NF;i++){A[$i]=i;}}else{a=A["20002:1465"];b=A["20004:1318"];c=A["20004:1319"];if ($a=="1" && ($b=="1" || $c=="1")){print $1;}}}'| gzip - -c >"$outfile" 
	;;
esac

exit 0
