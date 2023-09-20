#!/usr/bin/python3

import argparse
import pandas as pd
import numpy as np
import csv
import sys
import datetime
from functools import reduce

parser=argparse.ArgumentParser()
parser.add_argument('-i','--input',required=True,action='append',help="Input files")
parser.add_argument('-o','--output',required=True,action='store',help="Output prefix")

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    sys.exit(1)

infiles=args.input
out_prefix=args.output

print("input: %s\noutput prefix: %s\n" %(",".join(infiles),out_prefix),file=sys.stderr)

if len(infiles)==1:
    print("ERROR: only one input file provided",file=sys.stderr)
    sys.exit(1)

if len(infiles)>2:
    print("WARN: more than two input files provided; only merging the first two",file=sys.stderr)

inputDF=list()
for f in infiles:
    inputDF.append(pd.read_table(f,sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False))
    if len(inputDF)==2:
        break

print("Checking IDs in input DFs\n",file=sys.stderr)
for i in range(len(inputDF)):
    print("%s rows: %d" % (infiles[i],len(inputDF[i])),file=sys.stderr)
    print("%s columns: %d\n" % (infiles[i],len(inputDF[i].columns.values.tolist())),file=sys.stderr)
print("Done\n",file=sys.stderr)

merged=pd.merge(inputDF[0],inputDF[1],on="f.eid",how="outer",indicator="indicator",suffixes=("_left","_right"))
merged.replace(np.nan,"NA",inplace=True)
merged.replace("","NA",inplace=True)
L=list()
Lkeep=list()
for c in merged.columns.values.tolist():
    if c.endswith("_left"):
        x=c[:-5]
        L.append(x)
        merged[x]="NA"
    elif not c.endswith("_right") and c!="indicator":
        Lkeep.append(c)

print("INFO: intersecting columns: %d" %(len(L)),file=sys.stderr)
for c in L:
    print("INFO: current column: %s" %(c),file=sys.stderr)
    flag=True
    for i,row in merged.iterrows():
        if merged.at[i,"indicator"]=="both":
            if merged.at[i,c+"_left"]==merged.at[i,c+"_right"]:
                merged.at[i,c]=merged.at[i,c+"_left"]
            else:
                print("ERROR: %s" %(c))
                flag=False
                break
        elif merged.at[i,"indicator"]=="left_only":
            merged.at[i,c]=merged.at[i,c+"_left"]
        else:
            merged.at[i,c]=merged.at[i,c+"_right"]
    if flag:
        Lkeep.append(c)

print("INFO: writing output",file=sys.stderr)
merged.to_csv(out_prefix+".txt.gz",sep="\t",index=False,quotechar='"',quoting=csv.QUOTE_NONE,columns=Lkeep)

sys.exit(0)
