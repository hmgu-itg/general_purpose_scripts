#!/usr/bin/env bash

infile=$1
outfile=$2

tfile=${outfile/%.gz}

bcftools view -h "$infile" -Ov > "$tfile"
bcftools view -H "$infile" -Ov | perl -lne 'BEGIN{$n=0;}{$x=rand();if ($x<0.1){print $_;$n++;} exit 0 if $n>1000;}' >> "$tfile"
bgzip "$tfile"
tabix -p vcf "$outfile"
