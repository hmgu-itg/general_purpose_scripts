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
    echo "Usage: hesin_release.sh -i <main HESIN table> -d <HESIN DIAG table> -p <HESIN OPER table> -r <release> -o <output dir> { -x <list of sample IDs to exclude> -k <keep temp files>}"
    echo ""
    echo "All files are tab-separated."
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

outdir=""
release=""
exfile=""
keep="NO"
while getopts "hi:d:x:p:r:o:k" opt; do
    case $opt in
        i)main_fname=($OPTARG);;
        d)diag_fname=($OPTARG);;
        x)exfile=($OPTARG);;
        p)oper_fname=($OPTARG);;
        r)release=($OPTARG);;
        o)outdir=($OPTARG);;
        k)keep="YES";;
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
tmpdir=$(mktemp -d -p "$outdir" tempdir_hesin_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary directory in $outdir"
    exit 1
fi

: > $logfile

#-----------------------------------------------------------------------------------------------------

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
echo "TEMP DIR: $tmpdir" | tee -a "$logfile"
echo "KEEP TEMP FILES: $keep" | tee -a "$logfile"
echo "" | tee -a "$logfile"

#-----------------------------------------------------------------------------------------------------

datestr=$(date +%d-%b-%Y)
echo $release > "$tmpdir"/RELEASE
echo $datestr > "$tmpdir"/CREATED

#-----------------------------------------------------------------------------------------------------

declare -a ar
ar=("hesin.txt" "hesin_diag.txt" "hesin_oper.txt")
for fname in "$main_fname" "$diag_fname" "$oper_fname";do
    echo "REPLACING EMPTY FIELDS WITH NAs IN $fname " | tee -a "$logfile"
    outf="$tmpdir"/"${ar[0]}"
    echo "SAVING RESULT IN" "$outf" | tee -a "$logfile"
    echo "" | tee -a "$logfile"
    cat "$fname" | "$replace_script" > "$outf"
    ar=("${ar[@]:1}")
done

#-----------------------------------------------------------------------------------------------------

for fname in hesin.txt hesin_diag.txt hesin_oper.txt;do
    echo "EXCLUDING SAMPLES FROM $fname" | tee -a "$logfile"
    f="$tmpdir"/"$fname"
    eid_col=$(getColNum "$f" "eid" "cat")
    cat "$f" | perl -slne 'BEGIN{%h=();if (length($f)!=0){open(fh,"<",$f);while(<fh>){chomp;$h{$_}=1;}close(fh);}}{@a=split(/\t/);if (!defined($h{$a[$c-1]})){print $_;}}' -- -c="$eid_col" -f="$exfile" | sponge "$f"
done

#-----------------------------------------------------------------------------------------------------

echo "" | tee -a "$logfile"
echo "CREATING OUTPUT FILE" | tee -a "$logfile"
cd "$tmpdir" && tar -zcf "$outfile" hesin.txt hesin_diag.txt hesin_oper.txt RELEASE CREATED && cd -

#-----------------------------------------------------------------------------------------------------

if [[ "$keep" == "NO" ]];then
    rm -rf "$tmpdir"
fi

echo "" | tee -a "$logfile"
date "+%F %H-%M-%S"|tee -a "$logfile"
exit 0
