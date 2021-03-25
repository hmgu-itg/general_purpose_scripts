#!/usr/bin/env bash

indir=$1
threads=$2
c=$3

indir=${indir%/}

flist=$(mktemp -t merge_list_chr"$c"-XXXXXXX)
if [[ ! -f $flist ]];then
    echo "ERROR: could not create list file; exit"
    exit 1
fi

suffix=""
for f in $(find "$indir" -maxdepth 1 -mindepth 1 -name "*.vcf.gz" ! -name "merged*.vcf.gz");do
    realpath $f >> "$flist"
    if [[ -z "$suffix" ]];then
	b=$(basename $f)
	b=${b/#*filtered.}
	suffix=_chr"$c"_"$b"
    fi
done

output="$indir"/merged"$suffix".vcf.gz
plout="$indir"/merged"$suffix"

echo "MERGING VCF FILES AND CONVERTING TO PLINK FORMAT"
echo "INPUT DIR $indir"
echo "USING THREADS $threads"
echo "OUTPUT VCF $output"
echo "PLINK OUTPUT PREFIX $plout"
echo "START" $(date)
echo "-------------------------------------"
echo ""

# bcftools merge
if [[ ! -s "$output" ]];then
    echo $(date) "MERGING FILES:"
    cat "$flist"
    echo ""
    singularity exec -B /compute/Genomics /compute/Genomics/containers/worker_3.1 bcftools concat -a -f "$flist" -Oz -o "$output" --threads "$threads"
    if [[ $? -ne 0 ]];then
	echo "ERROR: merging failed; exit"
	if [[ -f "$output" ]];then rm -v "$output";fi
	rm -v "$flist"
	exit 1
    fi
    echo "MERGING DONE" $(date)
    echo ""
else
    echo "INFO: merged VCF $output already exists; skipping merging"
fi

# tabix
if [[ ! -s "$output".tbi ]];then
    singularity exec -B /compute/Genomics /compute/Genomics/containers/worker_3.1 tabix "$output"
    if [[ $? -ne 0 ]];then
	echo "ERROR: tabix failed; exit"
	if [[ -f "$output" ]];then rm -v "$output";fi
	if [[ -f "$output".tbi ]];then rm -v "$output".tbi;fi
	rm -v "$flist"
	exit 1
    fi
    echo "INDEXING DONE" $(date)
    echo ""
else
    echo "INFO: index file ${output}.tbi already exists; skipping indexing"
fi

# PLINK
flist2=$(mktemp -t plink_list_chr"$c"-XXXXXXX)
if [[ ! -f $flist2 ]];then
    echo "ERROR: could not create PLINK2 merging list file; exit"
    exit 1
fi

echo "CONVERTING CHUNKS TO PGEN"
cat "$flist" | while read fname
do
    outname=${fname/%.vcf.gz}
    /compute/Genomics/software/plink2/plink2_23.Mar.2021/plink2 --vcf "$fname" --make-pgen --out "$outname" --threads "$threads" --vcf-half-call m  
    echo "$outname".pgen "$outname".pvar "$outname".psam >> "$flist2"
done
echo "CONVERTING DONE"
echo ""

echo "MERGING PGEN CHUNKS"
/compute/Genomics/software/plink2/plink2_23.Mar.2021/plink2 --pmerge-list "$flist2" --make-pgen --out "$plout" --threads "$threads"
echo "DONE"
echo ""

# delete VCF chunk files
if [[ $? -eq 0 ]];then
    cat "$flist" | while read fname
    do
	rm -v "$fname"
    done
fi

# delete PGEN chunk files
cat "$flist" | while read fname
do
    outname=${fname/%.vcf.gz}
    rm -v "$outname".pgen "$outname".pvar "$outname".psam
done

rm -v "$flist"
rm -v "$flist2"

echo "INFO: done" $(date)
