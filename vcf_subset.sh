#!/usr/bin/env bash

infile=$1
outfile=$1

tfile=$(mktemp -t header-XXXXXXX)

bcftools view -h "$infile" -Ov > "$tfile"
bcftools view -H "$infile" -Ov | perl -lne 'BEGIN{$n=0;}{$x=rand();if ($x<0.1){print $_;$n++;} exit 0 if $n>1000;}' >> "$tfile"
bgzip "$tfile" > "$outfile"
tabix -p vcf "$outfile"
rm "$tfile"
