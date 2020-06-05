#!/usr/bin/bash

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

out=${out/%/}

mkdir -p "$out"

logfile="$out"/"cojo-wrapper.log"

echo "bfile         : $bfile" >> "$logfile"
echo "input table   : $input" >> "$logfile"
echo "output dir    : $out" >> "$logfile"

# reading the input table
cut -f 1 "$input"|sort|uniq| while read snp; do

echo "=========================================" >> "$logfile"
echo "Variant: $snp" >> "$logfile"

grep "^$snp" "$input"|while read id chr pos a1 a2 f1 b se p N

outKnown="$out"/"$suffix"."known.txt"

# sort if there are more than one ID in the uniprot string
uniprot=$(echo "$uniprot"| perl -lne '@a=split(/,/);print join(",",sort @a);')

start=$((pos-window))
end=$((pos+window))
intersectBed -wb -a <(echo "\"$chr $start $end\""| tr ' ' '\t') -b "$known" | awk -v x="$uniprot" 'BEGIN{FS="\t";OFS="\t";}$8==x{print $0;}' | cut -f 1,3| tr '\t' ':' > "$outKnown"

# if there are no known signals in the bp window
nKnown=$(cat "$outKnown"| wc -l)
if [[ "nKnown" -eq 0 ]];then
    echo "$id : no known signals" >> "$logfile"
    continue
fi

# if the tested variant is known
c=$(grep -c "$id" "$outKnown")
if [[ "$c" -ne 0 ]];then
    echo "$id : is known" >> "$logfile"
    continue
fi

# extract variants
plinkout="$out"/"$suffix"
echo -n "Extracting $id and known signals ... " >> "$logfile"
plink --make-bed --bfile "$bfile" --out "$plinkout" --extract <(echo "$id"|cat - "$outKnown") --allow-no-sex
echo "Done " >> "$logfile"

# if there are enough variants
c=$(cat "$plinkout".bim | wc -l)
if [[ "$c" -lt 2 ]];then
    echo "$id : not enough variants" >> "$logfile"
    continue
fi

# our variant should be in the PLINK output
c=$(cat "$plinkout".bim | grep "$id" | wc -l)
if [[ "$c" -eq 0 ]];then
    echo "$id : ERROR : not in bfile" >> "$logfile"
    continue
fi

# known signals which are not in the bfile
echo -n "$id: known signals not in the bfile: " >> "$logfile"
cut -f 2 "$plinkout"."bim" | cat - "$outKnown" | sort|uniq -u  >> "$logfile"

# create cojo file with m/a results
cojofile="$out"/"$suffix"."ma"
echo "$id $a1 $a2 $f1 $b $se $p $N" | tr ' ' '\t' > "$cojofile"
cat "$outKnown"| tr ':' ' '|while read cr ps;do
    tabix "$ma_path/$panel/$prot/$panel.$prot.metal.bgz" $cr:$ps-$ps| cut -f 1-5,9-11,17| awk -v c=$cr -v p=$ps 'BEGIN{FS="\t";OFS="\t";}$1==c && $2==p{print $1":"$2,$3,$4,$5,$6,$7,$8,$9;}'  >> "$cojofile"
done

# calling GCTA
echo -n "Calling GCTA ... " >> "$logfile"
gcta64 --bfile "$plinkout" --cojo-file "$cojofile" --cojo-cond "$outKnown" --out "$out"/"$suffix"."out"
echo "Done " >> "$logfile"

done

