# general_purpose_scripts
General purpose scripts
# cojo-prepare.sh

Given a table with meta-analysis results for independent hits, create an extended table 
containing in addition known signals in a specified bp window around hits.

`$ cojo-prepare.sh -i <input table> -o <output table> -k <known signals> -w <bp window>`
  
+ input table: a tab separated table containing meta-analysis results for independent hits; should contain following columns:
    +  column 1: panel
    +  column 2: protein
    +  column 5: uniprot ID of the protein
    +  column 6: chromosome of the hit
    +  column 8: position of the hit
    +  column 10: effect allele
    +  column 11: other allele
    +  column 12: frequency of the effect allele
    +  column 16: beta
    +  column 17: standard error
    +  column 18: p-value
    +  column 24: number of missing samples in meta-analysis of the current hit
+ output table: a tab separated table containing
    + record ID of the form: PANEL_PROTEIN_CHR_POS, where CHR and POS refer to the genomic position of the independent hit from the input table
    + CHR: chromosome of the signal (either input hit, or conditioning variant)
    + POS: position of the signal (either input hit, or conditioning variant)
    + A1
    + A2
    + freq
    + beta
    + SE
    + p-value
    + N
    + comma-separated list of samples used in the meta-analysis
+ known signals: a BED file containing information about known associations
+ optional bp window restricting conditioning known signals (default: 1Mbp)
