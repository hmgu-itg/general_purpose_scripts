#!/usr/bin/bash

function usage {
    echo ""
    echo "Usage: $0 -f <PLINK file prefix>"
    echo "          -o <output prefix>"
    echo "          -i <input table>"
    echo "          -w <bp window>"
    echo "          -k <BED file with known signals>"
    exit 0
}


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

