#!/usr/bin/python3

import argparse
import pandas as pd
import numpy as np
import csv
import sys
import datetime
import logging
import gc

# if a sample occurs in both dataframes and the values for that sample are different,
# then set the resulting value to np.nan
def fill_column(df,colname):
    c1=colname+"_left"
    c2=colname+"_right"
    conditions=[(df[c1]==df[c2]) & (df["indicator"]=="both"),df["indicator"]=="left_only",df["indicator"]=="right_only"]
    choices=[df[c1],df[c1],df[c2]]
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
LOGGER.addHandler(logging.StreamHandler(sys.stderr))

LOGGER.info("input: %s\noutput prefix: %s\n" %(",".join(infiles),out_prefix))

if len(infiles)==1:
    LOGGER.error("only one input file provided")
    sys.exit(1)

if len(infiles)>2:
    LOGGER.warning("more than two input files provided; only checking the first two")

columns=list()
for f in infiles:
    df=pd.read_table(f,sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,nrows=1)
    columns.append(set([x for x in df.columns.values.tolist() if x !="f.eid"]))
    if len(columns)==2:
        break
    
I=columns[0].intersection(columns[1])
if len(I)==0:
    LOGGER.info("no common columns")
    print("-i %s -i %s" %(infiles[0],infiles[1]))
    sys.exit(0)

C1=columns[0].difference(columns[1])
C2=columns[1].difference(columns[0])

C1.add("f.eid")
C2.add("f.eid")
I.add("f.eid")

DF1=pd.read_table(infiles[0],sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=list(I))
DF2=pd.read_table(infiles[1],sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=list(I))
merged=pd.merge(DF1,DF2,on="f.eid",how="outer",indicator="indicator",suffixes=("_left","_right"))
merged.replace(np.nan,"NA",inplace=True)
merged.replace("","NA",inplace=True)
for c in I:
    if c=="f.eid":
        continue
    merged[c]="NA"
for c in I:
    if c=="f.eid":
        continue
    fill_column(merged,c)
    # LOGGER.info("filled %s" %(c))

Lkeep=list()
Lkeep.append("f.eid")
L=list()
Lfile=list()
for c in I:
    if c=="f.eid":
        continue
    i=merged[c].isna().sum()
    if i>0:
        LOGGER.warning("%d mismatches in %s" %(i,c))
        L.append(c)
    else:
        Lkeep.append(c)

if len(L)!=0:
    for c in L:
        LOGGER.error("mismatches in %s" %(c))
    sys.exit(1)

fname=out_prefix+".txt.gz"
LOGGER.info("writing %d common columns to %s" %(len(I)-1,fname))
merged.to_csv(fname,sep="\t",index=False,quotechar='"',quoting=csv.QUOTE_NONE,columns=Lkeep)
Lfile.append("-i "+fname)

del DF1
del DF2
del merged
gc.collect()

if len(C1)!=0:
    fname=out_prefix+"1.txt.gz"
    LOGGER.info("writing %d columns to %s" %(len(C1)-1,fname))
    DF1=pd.read_table(infiles[0],sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=list(C1))
    DF1.to_csv(fname,sep="\t",index=False,quotechar='"',quoting=csv.QUOTE_NONE)
    Lfile.append("-i "+fname)
    
if len(C2)!=0:
    fname=out_prefix+"2.txt.gz"
    LOGGER.info("writing %d columns to %s" %(len(C2)-1,fname))
    DF1=pd.read_table(infiles[1],sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=list(C2))
    DF1.to_csv(fname,sep="\t",index=False,quotechar='"',quoting=csv.QUOTE_NONE)
    Lfile.append("-i "+fname)

print("%s" %(" ".join(Lfile)))

sys.exit(0)
