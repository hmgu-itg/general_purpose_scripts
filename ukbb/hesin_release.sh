#!/usr/bin/env bash

# ALL INPUT FILES ARE TAB SEPARATED

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
replace_script="${scriptdir}/replaceEmpty.pl"

function usage () {
    echo ""
    echo "Create a release for inpatient UKBB data"
    echo ""
    echo "Usage: hesin_release.sh -i <main HESIN table> -d <HESIN DIAG table> -p <HESIN OPER table> -r <release> -x <list of sample IDs to exclude> -o <output dir>"
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
exfile=""
while getopts "hi:d:x:p:r:o:" opt; do
    case $opt in
        i)main_fname=($OPTARG);;
        d)diag_fname=($OPTARG);;
        x)exfile=($OPTARG);;
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

if [[ ! -z "$exfile" ]];then
    exitIfNotFile "$exfile" "ERROR: exlude list $exfile does not exist"
fi

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
echo "IDs TO EXCLUDE: $exfile" | tee -a "$logfile"
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
echo "REPLACING EMPTY FIELDS WITH NAs IN MAIN" | tee -a "$logfile"
cat "$main_fname" |parallel --block 100M --pipe -N100000 "$replace_script" > "$tmpdir"/hesin.txt
echo "REPLACING EMPTY FIELDS WITH NAs IN DIAG" | tee -a "$logfile"
cat "$diag_fname" |parallel --block 100M --pipe -N100000 "$replace_script" > "$tmpdir"/hesin_diag.txt
echo "REPLACING EMPTY FIELDS WITH NAs IN OPER" | tee -a "$logfile"
cat "$oper_fname" |parallel --block 100M --pipe -N100000 "$replace_script" > "$tmpdir"/hesin_oper.txt

if [[ ! -z "$exfile" ]];then
    echo "EXCLUDING SAMPLES FROM $exfile IN MAIN" | tee -a "$logfile"
    tmp_fname="$tmpdir"/hesin.txt
    eid_col=$(getColNum "$tmp_fname" "eid" "cat")
    ncols=$(head -n 1 "$tmp_fname"| tr '\t' '\n'| wc -l)
    fmtstr="2.1"
    for i in $(seq 1 $ncols);do
	fmtstr="${fmtstr}"",1.$i"
    done
    echo "EID COLUMN: $eid_col" | tee -a "$logfile"
    echo "TOTAL COLUMNS: $ncols" | tee -a "$logfile"
    echo "FORMAT STRING: $fmt_str" | tee -a "$logfile"
    cat <(head -n 1 "$tmp_fname") <(join -1 "$eid_col" -2 1 -e NULL -a 1 -o "$fmtstr" -t$'\t' <(tail -n +2 "$tmp_fname"| sort -k"$eid_col","$eid_col") <(sort "$exfile")| grep NULL| cut --complement -f 1) | sponge "$tmp_fname"
    
    echo "EXCLUDING SAMPLES FROM $exfile IN MAIN" | tee -a "$logfile"
    tmp_fname="$tmpdir"/hesin_diag.txt
    eid_col=$(getColNum "$tmp_fname" "eid" "cat")
    ncols=$(head -n 1 "$tmp_fname"| tr '\t' '\n'| wc -l)
    fmtstr="2.1"
    for i in $(seq 1 $ncols);do
	fmtstr="${fmtstr}"",1.$i"
    done
    echo "EID COLUMN: $eid_col" | tee -a "$logfile"
    echo "TOTAL COLUMNS: $ncols" | tee -a "$logfile"
    echo "FORMAT STRING: $fmt_str" | tee -a "$logfile"
    cat <(head -n 1 "$tmp_fname") <(join -1 "$eid_col" -2 1 -e NULL -a 1 -o "$fmtstr" -t$'\t' <(tail -n +2 "$tmp_fname"| sort -k"$eid_col","$eid_col") <(sort "$exfile")| grep NULL| cut --complement -f 1) | sponge "$tmp_fname"
    
    echo "EXCLUDING SAMPLES FROM $exfile IN MAIN" | tee -a "$logfile"
    tmp_fname="$tmpdir"/hesin_oper.txt
    eid_col=$(getColNum "$tmp_fname" "eid" "cat")
    ncols=$(head -n 1 "$tmp_fname"| tr '\t' '\n'| wc -l)
    fmtstr="2.1"
    for i in $(seq 1 $ncols);do
	fmtstr="${fmtstr}"",1.$i"
    done
    echo "EID COLUMN: $eid_col" | tee -a "$logfile"
    echo "TOTAL COLUMNS: $ncols" | tee -a "$logfile"
    echo "FORMAT STRING: $fmt_str" | tee -a "$logfile"
    cat <(head -n 1 "$tmp_fname") <(join -1 "$eid_col" -2 1 -e NULL -a 1 -o "$fmtstr" -t$'\t' <(tail -n +2 "$tmp_fname"| sort -k"$eid_col","$eid_col") <(sort "$exfile")| grep NULL| cut --complement -f 1) | sponge "$tmp_fname"
fi

cd "$tmpdir" && tar -zcf "$outfile" hesin.txt hesin_diag.txt hesin_oper.txt RELEASE CREATED && cd -

rm -rf "$tmpdir"
date "+%F %H-%M-%S"|tee -a "$logfile"
