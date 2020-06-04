#!/usr/bin/bash

function usage {
    echo ""
    echo "Usage: $0 -f <PLINK file prefix>"
    echo "          -o <output dir>"
    echo "          -i <input table>"
    echo "          -w <bp window>"
    echo "          -k <BED file with known signals>"
    exit 0
}

# default
window=1000000

OPTIND=1
while getopts "f:o:i:w:k:" optname; do
    case "$optname" in
        "f" ) bfile="${OPTARG}";;
        "o" ) out="${OPTARG}";;
        "i" ) input="${OPTARG}";;
        "w" ) window="${OPTARG}";;
        "k" ) known="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

# total samples
totalN=$(cat "$bfile"."fam"| wc -l)

out=${out/%/}

mkdir -p "$out"

logfile="$out"/"cojo-wrapper.log"

cut -f 5,6,8,10-12,16-18,24 "$input"| while read uniprot chr pos a1 a1 f1 b se p nMiss; do
id="$chr:$pos"
suffix="$chr"_"$pos"
N=$((totalN-nMiss))

outKnown="$out"/"$suffix"."known.txt"

# sort if there are more than one IDs in the uniprot string
uniprot=$(echo "$uniprot"| perl -lne '@a=split(/,/);print join(",",sort @a);')

start=$((pos-window))
end=$((pos+window))
intersectBed -wb -a <(echo "\"$chr $start $end\""| tr ' ' '\t') -b "$known" | awk -v x="$uniprot" 'BEGIN{FS="\t";OFS="\t";}$8==x{print $0;}' | cut -f 1,3| tr '\t' ':' > "$outKnown"

# if there are no known signals in the bp window
nKnown=$(cat "$outKnown"| wc -l)
if [[ "nKnown" -eq 0 ]];then
    continue
fi

# if the tested variant is known
c=$(grep -c "$id" "$outKnown")
if [[ "$c" -ne 0 ]];then
    continue
fi

plinkout="$out"/"$suffix"

plink --make-bed --bfile "$bfile" --out "$plinkout" --extract <(echo "$id"|cat - "$outKnown") --allow-no-sex

# if there are enough variants
c=$(cat "$plinkout".bim | wc -l)
if [[ "$c" -lt 2 ]];then
    continue
fi

# make cojo file
cojofile="$out"/"$suffix"."ma"
echo "$id $a1 $a2 $f1 $b $se $p $N" | tr ' ' '\t' > "$cojofile"
cat "$outKnown"| while read x;do
    echo "$x NA NA NA NA NA NA NA" | tr ' ' '\t' >> "$cojofile"
done

#gcta64 --bfile "$plinkout" --cojo-file "$cojofile" --cojo-cond "$outKnown" --out "$out"/"$suffix"."out"

done

