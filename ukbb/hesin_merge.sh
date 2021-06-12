#!/usr/bin/env bash

# ALL INPUT FILES ARE TAB SEPARATED

bold=$(tput bold)
underlined=$(tput smul)
normal=$(tput sgr0)

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

if [[ -z "$release" ]];then
    echo "ERROR: release not specified" 1>&2
    exit 1
fi

if [[ -z "$outdir" ]];then
    echo "ERROR: output dir not specified" 1>&2
    exit 1
fi

if [[ ! -d $outdir ]];then
    echo "ERROR: output dir ($outdir) does not exist" 1>&2
    exit 1
fi

outdir=${outdir%/}
outfile="${outdir}/hesin_r${release}.txt.gz"
logfile="${outdir}/hesin_r${release}.log"
if [[ -f $outfile ]];then
    echo "ERROR: output file $outfile already exists" 1>&2
    exit 1
fi

: > $logfile

date | tee -a $logfile
echo "MAIN TABLE: $main_fname" | tee -a "$logfile"
echo "DIAG TABLE: $diag_fname" | tee -a "$logfile"
echo "OPER TABLE: $oper_fname" | tee -a "$logfile"
echo "RELEASE: $release" | tee -a "$logfile"
echo "OUTPUT DIR: $outdir" | tee -a "$logfile"
echo "OUTPUT FILE: $outfile" | tee -a "$logfile"
echo "" | tee -a "$logfile"

checkFields $main_fname "$logfile"
checkFields $diag_fname "$logfile"
checkFields $oper_fname "$logfile"
echo "" | tee -a "$logfile"

main_eid_col=$(getColNum $main_fname "eid")
if [[ -z $main_eid_col ]];then
    echo "ERROR: could not find \"eid\" column in $main_fname" 1>&2
    exit 1
fi
main_idx_col=$(getColNum $main_fname "ins_index")
if [[ -z $main_idx_col ]];then
    echo "ERROR: could not find \"ins_index\" column in $main_fname" 1>&2
    exit 1
fi
diag_eid_col=$(getColNum $diag_fname "eid")
if [[ -z $diag_eid_col ]];then
    echo "ERROR: could not find \"eid\" column in $diag_fname" 1>&2
    exit 1
fi
diag_idx_col=$(getColNum $diag_fname "ins_index")
if [[ -z $diag_idx_col ]];then
    echo "ERROR: could not find \"ins_index\" column in $diag_fname" 1>&2
    exit 1
fi
oper_eid_col=$(getColNum $oper_fname "eid")
if [[ -z $oper_eid_col ]];then
    echo "ERROR: could not find \"eid\" column in $oper_fname" 1>&2
    exit 1
fi
oper_idx_col=$(getColNum $oper_fname "ins_index")
if [[ -z $oper_idx_col ]];then
    echo "ERROR: could not find \"ins_index\" column in $oper_fname" 1>&2
    exit 1
fi

echo "MAIN TABLE eid COLUMN: $main_eid_col" | tee -a $logfile
echo "MAIN TABLE ins_index COLUMN: $main_idx_col" | tee -a $logfile
echo "DIAG TABLE eid COLUMN: $diag_eid_col" | tee -a $logfile
echo "DIAG TABLE ins_index COLUMN: $diag_idx_col" | tee -a $logfile
echo "OPER TABLE eid COLUMN: $oper_eid_col" | tee -a $logfile
echo "OPER TABLE ins_index COLUMN: $oper_idx_col" | tee -a $logfile
echo "" | tee -a $logfile

main_cols=$(head -n 1 $main_fname|tr '\t' '\n'|wc -l)
diag_cols=$(head -n 1 $diag_fname|tr '\t' '\n'|wc -l)
oper_cols=$(head -n 1 $oper_fname|tr '\t' '\n'|wc -l)

echo "COLUMNS IN MAIN TABLE: $main_cols" | tee -a $logfile
echo "COLUMNS IN DIAG TABLE: $diag_cols" | tee -a $logfile
echo "COLUMNS IN OPER TABLE: $oper_cols" | tee -a $logfile
echo "" | tee -a $logfile

main_str="1.1"
for i in $(seq 4 $((main_cols+1)));do
    main_str=$main_str",1.$i"
done

diag_str="2.4"
for i in $(seq 5 $((diag_cols+1)));do
    diag_str=$diag_str",2.$i"
done

join_str="1.1"
for i in $(seq 2 $((main_cols+diag_cols-3)));do
    join_str=$join_str",1.$i"
done

oper_str="2.4"
for i in $(seq 5 $((oper_cols+1)));do
    oper_str=$oper_str",2.$i"
done

# echo "MAIN STR: $main_str" | tee -a $logfile
# echo "DIAG STR: $diag_str" | tee -a $logfile
# echo "JOIN STR: $join_str" | tee -a $logfile
# echo "OPER STR: $oper_str" | tee -a $logfile

datestr=$(date +%F)

join --header -t $'\t' -1 1 -2 1 -e NA -a 1 -a 2 -o ${main_str},${diag_str} <(cat <(head -n 1 $main_fname|gawk -v n=$main_eid_col -v m=$main_idx_col 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i ~ /^ *$/){$i="NA";}}print $n"."$m,$0;}') <(tail -n +2 $main_fname|gawk -v n=$main_eid_col -v m=$main_idx_col 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i ~ /^ *$/){$i="NA";}}print $n"."$m,$0;}'|sort -k1,1)) <(cat <(head -n 1 $diag_fname|sed 's/level/diag_level/'|sed 's/arr_index/diag_arr_index/'|gawk -v n=$diag_eid_col -v m=$diag_idx_col 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i ~ /^ *$/){$i="NA";}}print $n"."$m,$0;}') <(tail -n +2 $diag_fname|gawk -v n=$diag_eid_col -v m=$diag_idx_col 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i ~ /^ *$/){$i="NA";}}print $n"."$m,$0;}'|sort -k1,1))| join --header -t $'\t' -1 1 -2 1 -e NA -a 1 -a 2 -o ${join_str},${oper_str} - <(cat <(head -n 1 $oper_fname|sed 's/level/oper_level/'|sed 's/arr_index/oper_arr_index/'|gawk -v n=$oper_eid_col -v m=$oper_idx_col 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i ~ /^ *$/){$i="NA";}}print $n"."$m,$0;}') <(tail -n +2 $oper_fname|gawk -v n=$oper_eid_col -v m=$oper_idx_col 'BEGIN{FS=OFS="\t";}{for (i=1;i<=NF;i++){if ($i ~ /^ *$/){$i="NA";}}print $n"."$m,$0;}'|sort -k1,1))|sed 's/\./\t/'|gawk -v d=$datestr -v r=$release 'BEGIN{FS=OFS="\t";}{if (NR==1){print $0,"CREATED","RELEASE";}else{print $0,d,r;}}' | gzip - -c > "$outfile"
