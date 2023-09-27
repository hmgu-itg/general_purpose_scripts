#!/usr/bin/python3

import argparse
import pandas as pd
import numpy as np
import csv
import sys
import datetime
import logging

# if a sample occurs in both dataframes and the values for that sample are different (and both != "NA")
# then set the resulting value to np.nan
def fill_column(df,colname):
    c1=colname+"_left"
    c2=colname+"_right"
    conditions=[(df[c1]==df[c2]) & (df["indicator"]=="both"),(df[c1]!=df[c2]) & (df[c1]=="NA") & (df["indicator"]=="both"),(df[c1]!=df[c2]) & (df[c2]=="NA") & (df["indicator"]=="both"),df["indicator"]=="left_only",df["indicator"]=="right_only"]
    choices=[df[c1],df[c2],df[c1],df[c1],df[c2]]
    df[colname]=np.select(conditions,choices,np.nan)

parser=argparse.ArgumentParser()
parser.add_argument('-i','--input',required=True,action='append',help="Input files")
parser.add_argument('-o','--output',required=True,action='store',help="Output prefix")
parser.add_argument("--verbose", "-v", help="Optional: verbosity level; default: info", required=False,choices=("debug","info","warning","error"),default="info")

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    sys.exit(1)

infiles=args.input
out_prefix=args.output
verbosity=logging.INFO
if args.verbose=="debug":
    verbosity=logging.DEBUG
elif args.verbose=="warning":
    verbosity=logging.WARNING
elif args.verbose=="error":
    verbosity=logging.ERROR

logfile=out_prefix+".log"
LOGGER=logging.getLogger("merge2")
LOGGER.setLevel(verbosity)
ch=logging.FileHandler(logfile,'w')
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(funcName)s -%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)
LOGGER.addHandler(logging.StreamHandler(sys.stdout))

LOGGER.info("input: %s\noutput prefix: %s\n" %(",".join(infiles),out_prefix))
#print("input: %s\noutput prefix: %s\n" %(",".join(infiles),out_prefix),file=sys.stderr)

if len(infiles)==1:
    LOGGER.error("only one input file provided")
    sys.exit(1)

if len(infiles)>2:
    LOGGER.warning("more than two input files provided; only merging the first two")

inputDF=list()
for f in infiles:
    inputDF.append(pd.read_table(f,sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False))
    if len(inputDF)==2:
        break

LOGGER.info("Checking IDs in input DFs\n")
for i in range(len(inputDF)):
    LOGGER.info("%s rows: %d" % (infiles[i],len(inputDF[i])))
    LOGGER.info("%s columns: %d\n" % (infiles[i],len(inputDF[i].columns.values.tolist())))
LOGGER.info("Done\n")

merged=pd.merge(inputDF[0],inputDF[1],on="f.eid",how="outer",indicator="indicator",suffixes=("_left","_right"))
merged.replace(np.nan,"NA",inplace=True)
merged.replace("","NA",inplace=True)

n_both=(merged.indicator.values=="both").sum()
n_left=(merged.indicator.values=="left_only").sum()
n_right=(merged.indicator.values=="right_only").sum()

LOGGER.info("total samples: %d" %(len(merged)))
LOGGER.info("common samples: %d" %(n_both))
LOGGER.info("samples in file 1 only: %d" %(n_left))
LOGGER.info("samples in file 2 only: %d" %(n_right))

L=list()
Lkeep=list()
for c in merged.columns.values.tolist():
    if c.endswith("_left"):
        x=c[:-5]
        L.append(x)
        merged[x]="NA"
    elif not c.endswith("_right") and c!="indicator":
        Lkeep.append(c)

LOGGER.info("intersecting columns: %d" %(len(L)))

for c in L:
    fill_column(merged,c)
    LOGGER.info("filled %s" %(c))
    
for c in L:
    i=merged[c].isna().sum()
    if i>0:
        LOGGER.warning("excluding %s (%d mismatches)" %(c,i))
    else:
        Lkeep.append(c)

LOGGER.info("writing output")
#merged.to_csv(out_prefix+".txt.gz",sep="\t",index=False,quotechar='"',quoting=csv.QUOTE_NONE,columns=Lkeep)
merged.to_csv(out_prefix+".txt.gz",sep="\t",index=False,quotechar='"',quoting=csv.QUOTE_NONE,na_rep="XXX")

sys.exit(0)
