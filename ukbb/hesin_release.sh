#!/usr/bin/env bash

# ALL INPUT FILES ARE TAB SEPARATED

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

function usage () {
    echo ""
    echo "Merging script for inpatient UKBB data"
    echo ""
    echo "Usage: hesin_merge.sh -i <main HESIN table> -d <HESIN DIAG table> -p <HESIN OPER table> -r <release> -o <output dir>"
    echo ""
    echo "Input files and output file are tab-separated."
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

outdir=""
release=""
while getopts "hi:d:p:r:o:" opt; do
    case $opt in
        i)main_fname=($OPTARG);;
        d)diag_fname=($OPTARG);;
        p)oper_fname=($OPTARG);;
        r)release=($OPTARG);;
        o)outdir=($OPTARG);;
        h)usage;;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

exitIfEmpty "$main_fname" "ERROR: main table (-i) not specified"
exitIfEmpty "$diag_fname" "ERROR: diag table (-d) not specified"
exitIfEmpty "$oper_fname" "ERROR: oper table (-p) not specified"
exitIfNotFile "$main_fname" "ERROR: main table $main_table does not exist"
exitIfNotFile "$diag_fname" "ERROR: diag table $diag_table does not exist"
exitIfNotFile "$oper_fname" "ERROR: oper table $oper_table does not exist"
exitIfEmpty "$release" "ERROR: release (-r) not specified"
exitIfEmpty "$outdir" "ERROR: output directory (-o) not specified"
exitIfNotDir "$outdir" "ERROR: $outdir is not a directory"

outdir=${outdir%/}
outfile="${outdir}/hesin_r${release}.tar.gz"
exitIfExists "$outfile" "ERROR: output file $outfile already exists"
logfile="${outdir}/hesin_r${release}.log"
outfile=$(realpath "$outfile")
logfile=$(realpath "$logfile")

: > $logfile

date "+%F %H-%M-%S"|tee -a "$logfile"
echo "Current dir: ${PWD}"|tee -a "$logfile"
echo "Command line: $scriptname ${args[@]}"|tee -a "$logfile"
echo "Output release: $release"|tee -a "$logfile"
echo ""|tee -a "$logfile"
echo "MAIN TABLE: $main_fname" | tee -a "$logfile"
echo "DIAG TABLE: $diag_fname" | tee -a "$logfile"
echo "OPER TABLE: $oper_fname" | tee -a "$logfile"
echo "OUTPUT DIR: $outdir" | tee -a "$logfile"
echo "OUTPUT FILE: $outfile" | tee -a "$logfile"
echo "" | tee -a "$logfile"

tmpdir=$(mktemp -d -p "$outdir" tempdir_hesin_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary directory in $outdir "|tee -a "$logfile"
    exit 1
fi

datestr=$(date +%d-%b-%Y)
echo $release > "$tmpdir"/RELEASE
echo $datestr > "$tmpdir"/CREATED

# replace empty fields with NAs
cat "$main_fname" | gawk -v 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i==""){$i="NA";}}print $n"."$m,$0;}' > "$tmpdir"/hesin.txt
cat "$diag_fname" | gawk -v 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i==""){$i="NA";}}print $n"."$m,$0;}' > "$tmpdir"/hesin_diag.txt
cat "$oper_fname" | gawk -v 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i==""){$i="NA";}}print $n"."$m,$0;}' > "$tmpdir"/hesin_oper.txt

cd "$tmpdir" && tar -zcf "$outfile" hesin.txt hesin_diag.txt hesin_oper.txt RELEASE CREATED && cd -

rm -rf "$tmpdir"
date "+%F %H-%M-%S"|tee -a "$logfile"
