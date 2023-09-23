#!/usr/bin/env bash

# ALL FILES ARE TAB SEPARATED

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

function usage () {
    echo ""
    echo "Add CREATED and RELEASE columns to input file"
    echo "Optionally exclude rows from input based on provided ID list"
    echo ""
    echo "Usage: add_release.sh -i <input.filename> -r <release> -o <output.filename> { -x <list of individual IDs to exclude> }"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

id_field="f.eid"
datestr=$(date +%d-%b-%Y)
infile=""
outfile=""
release=""
exclude_list=""

while getopts "hi:o:r:x:" opt; do
    case $opt in
        i)infile=($OPTARG);;
        o)outfile=($OPTARG);;
        r)release=($OPTARG);;
        x)exclude_list=($OPTARG);;
        h)usage;;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

exitIfEmpty "$infile" "ERROR: no input specified"
exitIfEmpty "$outfile" "ERROR: no output specified"
exitIfEmpty "$release" "ERROR: no release specified"
exitIfExists "$outfile" "ERROR: output file $outfile already exists"

# get zcat/cat command for the input file
cat=$(getCatCmd "$infile")
input_ID_column=$(getColNum "${infile}" "${id_field}" "${cat}")
exitIfEmpty "${input_ID_column}" "ERROR: no \"$id_field\" found in ${infile}"

#----------------------------------------------------

paste <("${cats}" "${infile}" | head -n 1) <(echo RELEASE CREATED | tr ' ' '\t') | gzip - -c > "${outfile}"
"${cat}" "${infile}" | tail -n +2 | perl -snle 'BEGIN{$,="\t";%h=();if (length($f)!=0){open(fh,"<",$f);while(<fh>){chomp;$h{$_}=1;}close(fh);}}{@a=split(/\t/);if (!defined($h{$a[$c-1]})){print $_,$r,$d;}}' -- -f="${exclude_list}" -c="${input_ID_column}" -r="${release}" -d="${datestr}" | gzip - -c >> "${outfile}"
    
exit 0

