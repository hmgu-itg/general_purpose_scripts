#!/bin/bash

function usage {
    echo ""
    echo "Usage: $0 -f <PLINK file prefix>"
    echo "          -o <output dir>"
    echo "          -i <input table>"
    exit 0
}

OPTIND=1
while getopts "f:o:i:" optname; do
    case "$optname" in
        "f" ) bfile="${OPTARG}";;
        "o" ) out="${OPTARG}";;
        "i" ) input="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

# required fields in the input table (tab separated)
# variant/protein association ID: panel_prot_chr_pos; chr and pos correspond to the variant being tested
#
# chr, pos, a1, a2, freq1, beta, SE, p, N: fields correspond to either the variant being tested, 
# or to conditioning variants, they have the same meaning as fields in a cojo file

out=${out/%/}

mkdir -p "$out"

logfile="$out"/"cojo-wrapper.log"

> "$logfile"

date >> "$logfile"
echo "bfile         : $bfile" >> "$logfile"
echo "input table   : $input" >> "$logfile"
echo "output dir    : $out" >> "$logfile"

# reading the input table
cut -f 1 "$input"|sort|uniq| while read varid; do

echo "=========================================" >> "$logfile"
echo "$varid" | tr '_' ' ' >> "$logfile"

id=$(echo $varid | cut -f 3- -d '_'| tr '_' ':')

# extract variants
plinkout="$out"/"$varid"
echo -n "Extracting $id and known signals ... " >> "$logfile"
plink --make-bed --bfile "$bfile" --out "$plinkout" --extract <(fgrep -w "$varid" "$input"|cut -f 2,3| tr '\t' ':') --allow-no-sex
echo "Done " >> "$logfile"
echo >> "$logfile"

# check if there are enough variants
c=$(cat "$plinkout".bim | wc -l)
if [[ "$c" -lt 2 ]];then
    echo "$varid : not enough variants" >> "$logfile"
    echo >> "$logfile"
    continue
fi

# our variant should be in the PLINK output
c=$(cat "$plinkout".bim | grep "$id" | wc -l)
if [[ "$c" -eq 0 ]];then
    echo "$id : ERROR : not in bfile" >> "$logfile"
    echo >> "$logfile"
    continue
fi

# known signals which are not in the bfile
echo "$id: known signals not in the bfile: " >> "$logfile"
cut -f 2 "$plinkout"."bim" | cat - <(fgrep -w "$varid" "$input"|cut -f 2,3| tr '\t' ':') | sort|uniq -u  >> "$logfile"
echo >> "$logfile"

# calling GCTA
echo -n "Calling GCTA ... " >> "$logfile"

echo "\nCOJO file: " >> "$logfile"
fgrep -w "$varid" "$input"|cut -f 2-|sed 's/\t/:/'| grep -v -f <(cut -f 2 "$plinkout"."bim" | cat - <(fgrep -w "$varid" "$input"|cut -f 2,3| tr '\t' ':') | sort|uniq -u) >> "$logfile"
echo "COND file" >> "$logfile"
fgrep -v -w "$id" "$plinkout"."bim" >> "$logfile"

gcta64 --bfile "$plinkout" --cojo-file <(fgrep -w "$varid" "$input"|cut -f 2-|sed 's/\t/:/'| grep -v -f <(cut -f 2 "$plinkout"."bim" | cat - <(fgrep -w "$varid" "$input"|cut -f 2,3| tr '\t' ':') | sort|uniq -u)) --cojo-cond <(fgrep -v -w "$id" "$plinkout"."bim") --out "$out"/"$varid"."out" 2>> "$logfile"
echo "Done " >> "$logfile"
echo >> "$logfile"

done
date >> "$logfile"

