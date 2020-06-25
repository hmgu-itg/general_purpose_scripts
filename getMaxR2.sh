#!/bin/bash

# get max r2 between the given variant and a list of valiants

function usage {
    echo ""
    echo "Usage: $0 -i <variant ID>"
    echo "          -l <list of variant IDs>"
    echo "          -f <pfile prefix>"
    exit 0
}


OPTIND=1
while getopts "i:l:f:" optname; do
    case "$optname" in
        "i" ) id="${OPTARG}";;
        "l" ) vars="${OPTARG}";;
        "f" ) pfile="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

ped="$pfile".ped
map="$pfile".map

n=$(fgrep -w -n $id $map|cut -d ':' -f 1)



# m/a results
ma_results="/storage/hmgu/projects/helic/OLINK/meta_analysis"

logfile="cojo-slct.log"

echo "" >> "$logfile"
date >> "$logfile"
echo "" >> "$logfile"

echo "ID            : $id" >> "$logfile"
echo "N             : $N" >> "$logfile"
echo "bfile         : $bfile" >> "$logfile"
echo "collinearity  : $ct" >> "$logfile"
echo "threads       : $threads" >> "$logfile"

read -r panel prot chr pos <<<$(echo $id|tr '_' ' ')

gcta64 --bfile "$bfile" --cojo-slct --out "$id".slct --cojo-file <(zcat "$ma_results"/"$panel"/METAL/"$panel"."$prot".metal.bgz| cut -f 1-5,9-11| awk -v n=$N 'BEGIN{FS="\t";OFS="\t";}{if (NR==1){print "SNP","A1","A2","freq","b","se","p","N";}else{ print $1":"$2,$3,$4,$5,$6,$7,$8,n;}}') --extract "$id".cond --threads "$threads" --cojo-collinear "$ct" --cojo-p "$pt"
