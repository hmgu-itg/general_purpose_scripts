#!/usr/bin/env bash

function usage () {
    echo "Convert UKBB enc to tab"
    echo ""
    exit 0
}

scriptname=$0
scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

config=""
infile=""
keyfile=""
outdir=""
execdir=""
encoding=""
while getopts "hk:o:e:c:i:" opt; do
    case $opt in
        i)infile=($OPTARG);;
        e)encoding=($OPTARG);;
        c)config=($OPTARG);;
        k)keyfile=($OPTARG);;
        o)outdir=($OPTARG);;
        h)usage;;
        *)usage;;
    esac
done

if [[ -z "$config" ]];then
    config="${scriptdir}/config.txt"
fi
exitIfEmpty "$outdir" "ERROR: output dir (-o) not specified"
exitIfEmpty "$keyfile" "ERROR: key file (-k) not specified"
exitIfEmpty "$infile" "ERROR: input (-i) not specified"
exitIfNotDir "$outdir" "ERROR: output dir $outdir is not a directory"
exitIfNotFile "$config" "ERROR: config file $config does not exist"
readValue "$config" EXEC_DIR execdir
execdir=${execdir%/}
exitIfNotDir "$execdir" "ERROR: EXEC_DIR $execdir is not a directory"
unpack_exe="$execdir"/ukbunpack
exitIfNotFile "$unpack_exe" "ERROR: $unpack_exe does not exist"
conv_exe="$execdir"/ukbconv
exitIfNotFile "$conv_exe" "ERROR: $conv_exe does not exist"
if [[ -z "$encoding" ]];then
    encoding="$execdir"/encoding.ukb
fi
exitIfNotFile "$encoding" "ERROR: encoding $encoding does not exist"

echo ""
echo "INFO: input: $infile"
echo "INFO: output dir: $outdir"
echo "INFO: key file: $keyfile"
echo "INFO: encoding: $encoding"
echo "INFO: config: $config"
echo ""

tmpd=$(mktemp -d -p "$outdir" temp_conv_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary dir in $outdir"
    exit 1
fi

cp "$encoding" "$tmpd" && cp "$unpack_exe" "$tmpd" && cp "$conv_exe" "$tmpd" && cp "$infile" "$tmpd" && cp "$keyfile" "$tmpd" && cd "$tmpd"
eval "./$(basename $unpack_exe) $(basename $infile) $(basename $keyfile)"
outfile1="$(basename $infile)"_"ukb"
exitIfNotFile "$outfile1" "ERROR: unpacked file $outfile1 does not exist"
eval "./$(basename $conv_exe) $outfile1 r"
outname=${outfile1/%enc_ukb/tab_dec}
rname=${outfile1/%enc_ukb/r}
Rscript <(cat "$rname" <(echo "write.table(bd,file=\"$outname\",quote=F,row.names=F,sep=\"\t\")"))
cp "$outname" "$outdir"

#rm -rf "$tmpd"

exit 0
