#!/usr/bin/env bash

infile=$1
outfile=$1

header=$(mktemp -t header-XXXXXXX)
body=$(mktemp -t body-XXXXXXX)

bcftools view -h "$infile" -Ov > "$header"
bcftools view -H "$infile" -Ov | perl -lne 'BEGIN{$n=0;}{$x=rand();if ($x<0.1){print $_;$n++;} exit 0 if $n>1000;}' > "$body"
cat "$header" "$body" | bgzip - > "$outfile"
tabix -p vcf "$outfile"
rm -v "$header" "$body"
