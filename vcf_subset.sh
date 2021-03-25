#!/usr/bin/env bash

infile=$1
outfile=$1

header=$(mktemp -t header-XXXXXXX)
body=$(mktemp -t body-XXXXXXX)

bcftools view -h "$infile" -Ov > "$header"
bcftools view -H "$infile" -Ov | shuf -n 1500 | sort -k2,2n > "$body"
cat "$header" "$body" | bgzip - > "$outfile"
tabix -p vcf "$outfile"
rm -v "$header" "$body"
