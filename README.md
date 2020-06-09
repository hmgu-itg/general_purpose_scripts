# general_purpose_scripts
General purpose scripts
# cojo-prepare.sh

Given a table with meta-analysis results for independent hits, create an extended table 
containing in addition known signals in a specified bp window around hits.

$ cojo-prepare.sh -i <input table> -o <output table> -k <known signals> -w <bp window: optional, default: 1000000>
  
+ input table: a tab separated table containing meta-analysis results for independent hits
+ output table: a tab separated table containing
    + record ID of the form: PANEL_PROTEIN_CHR_POS, where CHR and POS refer to the genomic position of the independent hit  
    from the input table
    + CHR: chromosome of the signal (either input hit, or conditioning variant)
    + POS: position of the signal (either input hit, or conditioning variant)
    + A1
    + A2
    + freq
    + beta
    + SE
    + p-value
    + N
    + comma-separated list of samples 
+ known signals
+ optional bp window restricting conditioning known signals (default: 1Mbp)
