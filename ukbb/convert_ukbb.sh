#!/usr/bin/env bash

function usage () {
    echo "Convert UKBB enc to tab file"
    echo ""
    echo "Usage: convert_ukbb.sh -i <input.enc>"
    echo "                       -k <input.key>"
    echo "                       -x <directory with UKBB executables>"
    echo "                       -o <output directory>"
    echo "                       -e <encoding.ukb; if not specified, assumed to be in \"-x\" directory>"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

scriptname=$0
scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

infile=""
keyfile=""
outdir=""
execdir=""
encoding=""
while getopts "hk:o:e:x:i:" opt; do
    case $opt in
        i)infile=($OPTARG);;
        e)encoding=($OPTARG);;
        x)execdir=($OPTARG);;
        k)keyfile=($OPTARG);;
        o)outdir=($OPTARG);;
        h)usage;;
        *)usage;;
    esac
done

exitIfEmpty "$outdir" "ERROR: output dir (-o) not specified"
exitIfEmpty "$keyfile" "ERROR: key file (-k) not specified"
exitIfEmpty "$infile" "ERROR: input (-i) not specified"
exitIfEmpty "$execdir" "ERROR: exec directory (-x) not specified"
exitIfNotDir "$execdir" "ERROR: exec directory $execdir is not a directory"
exitIfNotDir "$outdir" "ERROR: output dir $outdir is not a directory"
execdir=${execdir%/}
unpack_exe="$execdir"/ukbunpack
exitIfNotFile "$unpack_exe" "ERROR: $unpack_exe does not exist"
conv_exe="$execdir"/ukbconv
exitIfNotFile "$conv_exe" "ERROR: $conv_exe does not exist"
if [[ -z "$encoding" ]];then
    encoding="$execdir"/encoding.ukb
fi
exitIfNotFile "$encoding" "ERROR: encoding $encoding does not exist"
outdir=$(realpath "$outdir")

binput=$(basename "$infile")
outname=${binput/%enc/tab_dec}
exitIfExists "${outdir}/${outname}" "ERROR: output file $outname already exists"

echo ""
echo "INFO: input: $infile"
echo "INFO: output dir: $outdir"
echo "INFO: key file: $keyfile"
echo "INFO: encoding: $encoding"
echo "INFO: exec dir: $execdir"
echo ""

tmpd=$(mktemp -d -p "$outdir" temp_conv_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary dir in $outdir"
    exit 1
fi

{ cp "$encoding" "$tmpd" && cp "$unpack_exe" "$tmpd" && cp "$conv_exe" "$tmpd" && cp "$infile" "$tmpd" && cp "$keyfile" "$tmpd" && cd "$tmpd"; } || { echo "ERROR: could not copy files to $tmpd"; rm -rf "$tmpd"; exit 1; }
eval "./ukbunpack $(basename $infile) $(basename $keyfile)"
outfile1="$(basename $infile)"_"ukb"
exitIfNotFile "$outfile1" "ERROR: unpacked file $outfile1 does not exist"
eval "./ukbconv $outfile1 r"
rname=${outfile1/%enc_ukb/r}
tabname=${outfile1/%enc_ukb/tab}
exitIfNotFile "$rname" "ERROR: $rname does not exist"
exitIfNotFile "$tabname" "ERROR: $tabname does not exist"
Rscript <(cat "$rname" <(echo "write.table(bd,file=\"$outname\",quote=F,row.names=F,sep=\"\t\")"))
cp "$outname" "$outdir"

cd -
rm -rf "$tmpd"

exit 0
