#!/usr/bin/bash

function usage {
    echo ""
    echo "Usage: $0 -f <PLINK file prefix>"
    echo "          -o <output prefix>"
    echo "          -s <SNP to extract>"
    echo "          -c <file with conditioning SNPs>"
    echo "          -m <file with meta-analysis results>"
    exit 0
}


OPTIND=1
while getopts "f:o:s:c:m:" optname; do
    case "$optname" in
        "f" ) bfile="${OPTARG}";;
        "o" ) out="${OPTARG}";;
        "s" ) snp="${OPTARG}";;
        "c" ) cond="${OPTARG}";;
        "m" ) meta="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

gcta64 --bfile "$bfile" --out "$out" --extract-snp "$snp" --cojo-file "$meta" --cojo-cond "$cond"
