#!/usr/bin/python3

import argparse
import pandas as pd
import csv
import sys
import logging
import os
import re
import functools as ft

def transformExpr(string):
    return re.sub(r"(\w+)",lambda x:"(\'"+x.group(1)+"\' in L)" if (x.group(1).lower()!="and" and x.group(1).lower()!="or") else x.group(1).lower(),string)

def testF(L1,L2):
    return ft.reduce(lambda a,b:a or b,list(map(lambda x:eval(x,{},{"L":L1}),L2)))

verbosity=logging.INFO

parser=argparse.ArgumentParser()
parser.add_argument('--input','-i',required=True,action='store',help="Input file")
parser.add_argument('-icd10','--icd10',required=True,action='store',help="List of ICD10 codes")
parser.add_argument('-o','--output',required=True,action='store',help="Output file")
parser.add_argument("--verbose", "-v", help="Verbosity level; default: info",required=False,choices=("debug","info","warning","error"),default="info")

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    sys.exit(1)

infile=args.input
outfname=args.output

if args.verbose is not None:
    if args.verbose=="debug":
        verbosity=logging.DEBUG
    elif args.verbose=="warning":
        verbosity=logging.WARNING
    elif args.verbose=="error":
        verbosity=logging.ERROR

LOGGER=logging.getLogger("test_hesin_select")
LOGGER.setLevel(verbosity)
ch=logging.StreamHandler(sys.stderr)
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(funcName)s -%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

icd10codes=set()
if not args.icd10 is None:
    with open(args.icd10,"r") as f:
        for l in f.read().splitlines():
            if not re.match("^\s*$",l) and not re.match("^#.*",l):
                icd10codes.add(l)

LOGGER.info("input file: %s" % infile)
LOGGER.info("output file: %s" % outfname)

#-----------------------------------------------------------------------------------------------------------------------------

L=list()
df_diag=pd.read_table(infile,sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","ins_index","diag_icd9","diag_icd10"])
#print(df_diag)
df_diag2=df_diag.groupby(["eid","ins_index"],as_index=False).agg({"diag_icd9":lambda x:list(x),"diag_icd10":lambda x:list(x)})
#print(df_diag2)
#print(icd10codes)
#print([transformExpr(z) for z in icd10codes])
#df_diag2["new"]=df_diag2["diag_icd10"].apply(lambda x:testF(x,[transformExpr(z) for z in icd10codes]))
#df_diag2.to_csv("temp",sep="\t",index=False,quotechar='"',quoting=csv.QUOTE_NONE)
#print(df_diag2)

L=df_diag2[df_diag2.apply(lambda x:testF(x["diag_icd10"],[transformExpr(z) for z in icd10codes])==True,axis=1)]["eid"].unique().tolist()
LOGGER.info("output %d ID(s)" %len(L))
with open(outfname,"w") as f:
    if len(L):
        print("\n".join(L),file=f)
