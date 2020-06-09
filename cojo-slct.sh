#!/bin/bash

function usage {
    echo ""
    echo "Usage: $0 -i <signal ID>"
    echo "          -n <number of samples in m/a analysis>"
    echo "          -f <bfile prefix>"
    exit 0
}

OPTIND=1
while getopts "i:n:f:e:" optname; do
    case "$optname" in
        "i" ) id="${OPTARG}";;
        "n" ) N="${OPTARG}";;
        "f" ) bfile="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

# a simple wrapper for cojo-slct

# m/a results
ma_results="/storage/hmgu/projects/helic/OLINK/meta_analysis"

read -r panel prot chr pos <<<$(echo $id|tr '_' ' ')

gcta64 --bfile "$bfile" --cojo-slct --out "$id".slct --cojo-file <(zcat "$ma_results"/"$panel"/METAL/"$panel"."$prot".metal.bgz| cut -f 1-5,9-11| awk -v n=$N 'BEGIN{FS="\t";OFS="\t";}{if (NR==1){print "SNP","A1","A2","freq","b","se","p","N";}else{ print $1":"$2,$3,$4,$5,$6,$7,$8,n;}}') --extract "$id".cond
