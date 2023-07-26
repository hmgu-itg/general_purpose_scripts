#!/bin/bash

function usage {
    echo ""
    echo "Usage: $0 -s <SNP>"
    echo "          -i <SNP list>"
    echo "          -b <input PLINK prefix>"
    echo "          -o <output>"
    exit 0
}


OPTIND=1
while getopts "i:o:s:b:" optname; do
    case "$optname" in
        "i" ) input_list="${OPTARG}";;
        "o" ) output="${OPTARG}";;
        "s" ) input_snp="${OPTARG}";;
        "b" ) input_prefix="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

tmp_prefix=/tmp/tmp_extr
tmp_prefix2=/tmp/tmp_LD

plink --make-bed --bfile "${input_prefix}" --out "${tmp_prefix}" --extract <(cat "${input_list}" <(echo "${input_snp}")|sort|uniq)
n=$(cat "${input_list}" <(echo "${input_snp}")|sort|uniq| wc -l)
n=$(( n + 1 ))
start=$(cut -f 4 "${tmp_prefix}".bim|sort -n|head -n 1)
end=$(cut -f 4 "${tmp_prefix}".bim|sort -rn|head -n 1)
k=$(( end - start ))
k=$(( k / 1000 ))
k=$(( k + 1 ))

plink --bfile "${tmp_prefix}" --out "${tmp_prefix2}" --r2 --ld-window-r2 0  --ld-snp "${input_snp}" --ld-window $n --ld-window-kb $k
# sed 's/^  *//' "${tmp_prefix2}".ld|sed 's/  */ /g'|cut -d ' ' -f 6,7| tail -n +2|awk -v v="${input_snp}" '{if (v!=$1) {print $0;}}' > "${output}"
sed 's/^  *//' "${tmp_prefix2}".ld|sed 's/  */ /g'|cut -d ' ' -f 6,7| tail -n +2 > "${output}"

rm -f "${tmp_prefix}".bim "${tmp_prefix}".fam "${tmp_prefix}".bed "${tmp_prefix}".nosex "${tmp_prefix}".log "${tmp_prefix2}".log "${tmp_prefix2}".nosex "${tmp_prefix2}".ld

exit 0

