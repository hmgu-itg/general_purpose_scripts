#!/usr/bin/env bash

scriptname=$0
args=("$@")

declare -i rep=$1
outdir=$2

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
ginput="${scriptdir}/generate_input.sh"
gupdate="${scriptdir}/generate_update.sh"

tmpdir=$(mktemp -d -p "$outdir" tempdir_test_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary directory in $outdir"
    exit 1
fi







exit 0

