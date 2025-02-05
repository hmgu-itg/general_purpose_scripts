#!/usr/bin/env bash

# ALL FILES ARE TAB SEPARATED

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

function usage () {
    echo ""
    echo "Add/update CREATED and RELEASE information"
    echo "Optionally exclude rows from input based on provided ID list"
    echo ""
    echo "Usage: add_release.sh -i <input file>"
    echo "                    { -o <optional: output directory> }"
    echo "                    { -f <optional: ID field name; default: \"f.eid\"> }"
    echo "                    { -r <optional: output release> }"
    echo "                    { -b <optional: basename of the output file; default: \"phenotypes\"> }"
    echo "                    { -x <optional: list of individual IDs to exclude> }"
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
bname="phenotypes"
outdir=""

while getopts "hf:i:o:b:x:r:" opt; do
    case $opt in
        f)id_field=($OPTARG);;
        i)infile=($OPTARG);;
        o)outdir=($OPTARG);;
        b)bname=($OPTARG);;
        x)exclude_list=($OPTARG);;
        r)release=($OPTARG);;
        h)usage;;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

exitIfEmpty "$infile" "ERROR: no input specified"
if [[ -z "$outdir" ]];then
    outdir=$(dirname "$infile")
fi
exitIfEmpty "$outdir" "ERROR: no output dir specified"
exitIfNotDir "$outdir" "ERROR: output dir $outdir is not a directory"
if [[ ! -w "$outdir" ]];then
    echo "ERROR: output dir $outdir is not writable" 1>&2
    exit 1
fi
outdir=${outdir%/}

# get zcat/cat command for the input file
cat=$(getCatCmd "$infile")
input_ID_column=$(getColNum "${infile}" "${id_field}" "${cat}")
exitIfEmpty "${input_ID_column}" "ERROR: \"$id_field\" not found in ${infile}"

prev_release=$(getRelease "${infile}" "${cat}")
if [[ -z "$prev_release" ]];then
    if [[ -z "$release" ]];then
	echo "ERROR: could not find RELEASE in $infile; release (-r) not specified" 1>&2
	exit 1
    fi
else
    if [[ -z "$release" ]];then
	release=$((prev_release+1))
    else
	if [[ "$prev_release" == "$release" ]];then
	    echo "ERROR: previous release ($prev_release) and specified release ($release) are the same" 1>&2
	    exit 1
	fi
    fi
fi

outfile="${outdir}/${bname}_r${release}.txt.gz"
exitIfExists "$outfile" "ERROR: output file $outfile already exists"

logfile="${outdir}/${bname}_r${release}.log"
: > "$logfile"

echo "" | tee -a "$logfile"
date "+%F %H-%M-%S" | tee -a "$logfile"
echo "" | tee -a "$logfile"

echo -e "\n---------------------------------------------------------\n" | tee -a "$logfile"

echo "Current dir: ${PWD}" | tee -a "$logfile"
echo "Command line: $scriptname ${args[@]}" | tee -a "$logfile"
echo "" | tee -a "$logfile"
echo "INPUT FILE: ${infile}" | tee -a "$logfile"
echo "OUTPUT DIR: ${outdir}" | tee -a "$logfile"
echo "ID FIELD: $id_field" | tee -a "$logfile"
echo "OUTPUT RELEASE: $release" | tee -a "$logfile"
echo "OUTPUT FILE: $outfile" | tee -a "$logfile"
echo "EXCLUDE LIST: $exclude_list" | tee -a "$logfile"
echo "LOG FILE: $logfile" | tee -a "$logfile"

echo -e "\n---------------------------------------------------------\n" | tee -a "$logfile"

#----------------------------------------------------

if [[ -z "$prev_release" ]];then
    paste <("${cat}" "${infile}" | head -n 1) <(echo RELEASE CREATED | tr ' ' '\t') | gzip - -c > "${outfile}"
    "${cat}" "${infile}" | tail -n +2 | perl -snle 'BEGIN{$,="\t";%h=();if (length($f)!=0){open(fh,"<",$f);while(<fh>){chomp;$h{$_}=1;}close(fh);}}{@a=split(/\t/);if (!defined($h{$a[$c-1]})){print $_,$r,$d;}}' -- -f="${exclude_list}" -c="${input_ID_column}" -r="${release}" -d="${datestr}" | gzip - -c >> "${outfile}"
else
    rcol=$(getColNum "$infile" "RELEASE" "$cat")
    ccol=$(getColNum "$infile" "CREATED" "$cat")
    "${cat}" "${infile}" | head -n 1 | gzip - -c > "${outfile}"
    "${cat}" "${infile}" | tail -n +2 | perl -snle 'BEGIN{$,="\t";%h=();if (length($f)!=0){open(fh,"<",$f);while(<fh>){chomp;$h{$_}=1;}close(fh);}}{@a=split(/\t/);if (!defined($h{$a[$c-1]})){$a[$y-1]=$r;$a[$z-1]=$d;print join("\t",@a);}}' -- -y="$rcol" -z="$ccol" -f="${exclude_list}" -c="${input_ID_column}" -r="${release}" -d="${datestr}" | gzip - -c >> "${outfile}"
fi

date "+%F %H-%M-%S" | tee -a "$logfile"
    
exit 0

