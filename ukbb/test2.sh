#!/usr/bin/env bash

args=("$@")
scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

f1=${args[0]}
f2=${args[1]}
outfile=${args[2]}
chunks=${args[3]}

echo "f1: $f1"
echo "f2: $f2"
echo "output: $outfile"
echo "chunks: $chunks"

join_two_files "LOG" "/tmp" $f1 $f2 $outfile $chunks

exit 0
