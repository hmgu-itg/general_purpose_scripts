#!/usr/bin/bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $scriptname))
source "${scriptdir}/ukbb/functions.sh"

perl_exe=$(/usr/bin/which perl)
filter_script="${scriptdir}/filter_bim_input.pl"
dentist_exe="/lustre/groups/itg/shared/DENTIST.tmp2"

# ref panel prefix: /lustre/groups/itg/teams/data-informatics/projects/ukbb/genotypes/imputed/EGAD00010001474/GO2_var_ids/ukbb.imputed.v3

function usage () {
    echo ""
    echo "Running DENTIST on SUSIE input file"
    echo ""
    echo "Usage: dentist_run.sh -s <Signal.sumstats.input> -i <SUSIE.input.txt.gz> -o <output.dir> -p <LD.panel.prefix> -k <keep temp files> -t <threads; default: 1>"
    echo ""
    echo "SUSIE.input.txt.gz is output of the susie_prepare.sh script"
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

input_fname=""
sumstat_fname=""
output_dir=""
panel_prefix=""
keep="NO"
threads=1
while getopts "hi:s:o:t:p:k" opt; do
    case $opt in
        i)input_fname=($OPTARG);;
        s)sumstat_fname=($OPTARG);;
        o)output_dir=($OPTARG);;
        p)panel_prefix=($OPTARG);;
        t)threads=($OPTARG);;
	k)keep="YES";;
        h)usage;;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

exitIfEmpty "$input_fname" "ERROR: no input specified"
exitUnlessExists "$input_fname" "ERROR: $input_fname does not exist"
exitIfEmpty "$sumstat_fname" "ERROR: no input specified"
exitUnlessExists "$sumstat_fname" "ERROR: $input_fname does not exist"
exitIfEmpty "$output_dir" "ERROR: no output directory specified"
exitUnlessExists "$output_dir" "ERROR: $output_dir does not exist"
exitIfNotDir "$output_dir" "ERROR: $output_dir is not a directory"
exitIfEmpty "$panel_prefix" "ERROR: no LD panel prefix specified"

output_dir=${output_dir/%\/}

b=$(basename $input_fname)
out=${b/%.txt.gz}
gwas=${b/%.txt.gz/.gwas}

log_fname="${output_dir}"/"${out}".script.log

: > "${log_fname}"

echo "SUMSTAT: $sumstat_fname" | tee -a "${log_fname}"
echo "INPUT: $input_fname" | tee -a "${log_fname}"
echo "REF PANEL PREFIX: $panel_prefix" | tee -a "${log_fname}"
echo "OUTPUT DIR: $output_dir" | tee -a "${log_fname}"
echo "KEEP TEMP FILES: $keep" | tee -a "${log_fname}"
echo "THREADS: $threads" | tee -a "${log_fname}"
echo "PERL: $perl_exe" | tee -a "${log_fname}"

# extract variants
chr=$(echo $b|cut -d '_' -f 3| cut -d ':' -f 1)
echo "CHR: $chr" | tee -a "${log_fname}"
exitUnlessExists "${panel_prefix}.chr${chr}.bed" "ERROR: ${panel_prefix}.chr${chr}.bed does not exist"
exitUnlessExists "${panel_prefix}.chr${chr}.bim" "ERROR: ${panel_prefix}.chr${chr}.bim does not exist"
exitUnlessExists "${panel_prefix}.chr${chr}.fam" "ERROR: ${panel_prefix}.chr${chr}.fam does not exist"
temp_ex=$(mktemp -p "${output_dir}" dentist.ex.XXXXXX)
echo "TEMP FILE: $temp_ex" | tee -a "${log_fname}"
zcat "${input_fname}"|cut -d ' ' -f 1| tail -n +2 > "${temp_ex}"
n=$(cat "${temp_ex}"| wc -l)
echo "EXTRACTED $n VARIANTS TO ${temp_ex}" | tee -a "${log_fname}"
plink --bfile "${panel_prefix}".chr"${chr}" --out "${output_dir}"/"${out}" --extract "${temp_ex}" --make-bed >>"${log_fname}" 2>>"${log_fname}"

# creating GWAS input for DENTIST
echo "CREATING GWAS INPUT FOR DENTIST" | tee -a "${log_fname}"
"${perl_exe}" "${filter_script}" "${output_dir}"/"${out}".bim  "${sumstat_fname}"  > "${output_dir}"/"${gwas}"

# running DENTIST
dentist_out="${output_dir}"/"${out}".dentist
${dentist_exe} --gwas-summary "${output_dir}"/"${gwas}" --out "${dentist_out}" --bfile "${output_dir}"/"${out}" --thread-num "$threads" --wind-dist 1000000 --maf 0.00001 >>"${log_fname}" 2>>"${log_fname}"

# temp files
if [[ $keep == "NO" ]];then
    echo "DELETING TEMP FILES" | tee -a "${log_fname}"
    rm "${temp_ex}" "${output_dir}"/"${out}".log "${output_dir}"/"${out}".nosex "${output_dir}"/"${out}".bed "${output_dir}"/"${out}".bim "${output_dir}"/"${out}".fam  "${output_dir}"/"${gwas}"
fi

exit 0
