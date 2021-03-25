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
	suffix=_chr"$c"."$b"
	suffix=${suffix/%.vcf.gz}
    fi
done

output="$indir"/merged"$suffix".vcf.gz
plout="$indir"/merged"$suffix"

echo "MERGING VCF FILES AND CONVERTING TO PGEN FORMAT"
echo "INPUT DIR $indir"
echo "USING THREADS $threads"
echo "OUTPUT VCF $output"
echo "PLINK OUTPUT PREFIX $plout"
echo "START" $(date)
echo "-------------------------------------"
echo ""

newest_vcf=$(find "$indir" -maxdepth 1 -mindepth 1 -name "*.vcf.gz" ! -name "merged*.vcf.gz" -printf "%T@ %p\n"| sort -t ' ' -k1,1n |cut -d ' ' -f 2- | tail -n 1)

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
    if [[ "$output" -ot "$newest_vcf" ]];then
	echo "INFO: merged VCF $output is older than input VCFs; removing and re-merging"
	rm -v "$output"
	if [[ -f "$output".tbi ]];then rm -v "$output".tbi;fi
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
fi

# tabix
echo "INDEXING" $(date)
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
    if [[ "$output".tbi -ot "$output" ]];then
	echo "INFO: index file is older than merged VCF; re-indexing"
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
	rm "$fname" "$fname".tbi
    done
fi

# delete PGEN chunk files
cat "$flist" | while read fname
do
    outname=${fname/%.vcf.gz}
    rm "$outname".pgen "$outname".pvar "$outname".psam
done

rm "$indir"/*.log
rm "$plout"-merge.*

rm -v "$flist"
rm -v "$flist2"

echo "INFO: done" $(date)
