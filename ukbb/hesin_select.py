#!/usr/bin/python3

import argparse
import pandas as pd
import csv
import sys
import tarfile

parser=argparse.ArgumentParser()
parser.add_argument('-i','--input',required=True,action='store',help="Input file")
parser.add_argument('-icd9','--icd9',required=False,action='store',help="List of ICD9 codes")
parser.add_argument('-icd10','--icd10',required=False,action='store',help="List of ICD10 codes")
parser.add_argument('-o','--output',required=True,action='store',help="Output prefix")

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    sys.exit(1)

if args.icd9 is None and args.icd10 is None:
    print("ERROR: no ICD codes given",file=sys.stderr)
    sys.exit(1)
    
infile=args.input
out_prefix=args.output
logF=open(out_prefix+".log","w")

icd9codes=list()
icd10codes=list()

if not args.icd9 is None:
    with open(args.icd9,"r") as f:
        icd9codes=f.read().splitlines()

if not args.icd10 is None:
    with open(args.icd10,"r") as f:
        icd10codes=f.read().splitlines()

#-----------------------------------------------------------------------------------------------------------------------------

with tarfile.open(infile,"r:*") as tar:
    df=pd.read_table(tar.extractfile("hesin_diag.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","diag_icd9","diag_icd10"])
    if icd10codes:
        L=df.loc[df["diag_icd10"].isin(icd10codes)]["eid"].unique().tolist()
    if icd9codes:
        L.extend(x for x in df.loc[df["diag_icd9"].isin(icd9codes)]["eid"].unique().tolist() if not x in L)
    with open(out_prefix+".txt","w") as f:
        print(*L,sep="\n",file=f)
