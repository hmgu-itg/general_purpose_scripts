#!/usr/bin/Rscript

## Synchronises a phenotype file to the exact same samples as a VCF file

argv=commandArgs(T)
library(data.table)
library(Hmisc)

if(length(argv)<3){
  cat("ERROR: At least 3 arguments expected: vcf_file pheno_file output_file [optional_exclusion_file]\n")
  exit(1)
}

vcffile=argv[1]
if(!file.exists(vcffile)){
  cat("ERROR: the specified VCF does not exist.\n")
  exit(1)
}

phenofile=argv[2]
if(!file.exists(phenofile)){
  cat("ERROR: the specified phenotype file does not exist.\n")
  exit(1)
}

outfile=argv[3]

exclfile=NULL
if(length(argv)>3){
  exclfile=argv[4]
  if(!file.exists(phenofile)){
    cat("ERROR: Exclusion file specified but the file does not exist.\n")
    exit(1)
  }
}


## Get VCF samples
header=system(paste0("zgrep -m1 '#CHROM' ", vcffile), intern=T)
vcfsamples=strsplit(header, "\t")
vcfsamples=vcfsamples[[1]][10:length(vcfsamples[[1]])]

## Get exclusion
exclist=NULL
if(!is.null(exclfile)){
  exclist=fread(exfile, header=F)$V1
}

## read sample file
pheno=fread(phenofile)

### remove samples not in VCF
pheno=pheno[ID %in% c(vcfsamples, 0)]

### add samples in VCF but not in sample file
phenoname=colnames(pheno)[2]
setnames(pheno, phenoname, "pheno")
pheno=rbind(pheno, data.table(ID=vcfsamples[vcfsamples %nin% pheno$ID], pheno=NA))
setnames(pheno, "pheno", phenoname)

### synchronise
pheno=rbind(pheno[1], pheno[match(vcfsamples, pheno$ID)])

### subset if warranted
if(!is.null(exclist)){
  pheno=pheno[ID %nin% exclist]
}

fwrite(pheno, outfile, quote=F, sep="\t")
