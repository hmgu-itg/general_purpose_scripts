#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

function usage () {
    echo ""
    echo "Sample QC for UKBB directly typed genotypes"
    echo ""
    echo "Usage: ukbb_sampleQC.sh -i | --input <space separated file with bed,bim filenames per row>"
    echo "                        -f | --fam <input fam file>"
    echo "                        -o | --output <output dir>"
    echo "                        -q | --qc <file with sample qc statistics>"
    echo "                        -h | --help"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

OPTS=$(getopt -o hq:i:o: -l help,input:,qc:,output-dir: -n 'ukbb_sampleQC' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

outdir=""
infile=""
qcfile=""
famfile=""
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
    -i|--input ) infile=$2; shift 2 ;;
    -f|--fam ) famfile=$2; shift 2 ;;
    -o|--output-dir ) outdir=$2; shift 2 ;;
    -q|--qc ) qcfile=$2; shift 2 ;;
    --) shift ; break ;;
    * ) echo "Unexpected option: $1" ; usage ;;
  esac
done

#-------------------------------------------------------------------------------------------------------------
# file/dir checks

#-------------------------------------------------------------------------------------------------------------

if [[ $(cat "$famfile"| wc -l) -ne () ]];then
    echo "ERROR: $famfile and $qcfile have differemt number of rows"
    exit 1
fi

