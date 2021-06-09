#!/usr/bin/python3

import argparse
import pandas as pd
import csv
import sys
from functools import reduce


parser=argparse.ArgumentParser()
parser.add_argument('-i','--input',required=True,action='store',help="Input file")
parser.add_argument('-u','--update',required=True,action='append',help="Update file(s)")

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    sys.exit(0)

infile=args.input
updates=args.update

print("input: %s\nupdates: %s\n" %(infile,",".join(updates)))

inputDF=pd.read_table(infile,sep="\t",header=[0,1],dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False)
print(inputDF)
print(inputDF.columns.values.tolist())
for x in inputDF.columns.values.tolist():
    if not x[0] in ("CREATED","RELEASE","f.eid"):
        print("%s\t%s" % (x[0],x[1]))

updateDF=list()
for f in updates:
    updateDF.append(pd.read_table(f,sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False))

print("Checking IDs in update DFs")
for i in range(len(updateDF)):
    for j in range(len(updateDF)):
        if i<j:
            if not updateDF[i]["f.eid"].equals(updateDF[j]["f.eid"]):
                print("%d,%d not equal" %(i,j))

print("Checking columns in update DFs")
for i in range(len(updateDF)):
    for j in range(len(updateDF)):
        if i<j:
            if len([x for x in updateDF[i].columns.values.tolist() if x in updateDF[j].columns.values.tolist()])>1:
                print("%d,%d have common columns" %(i,j))
                
print("Merging update DFs")
merged=reduce(lambda x, y:pd.merge(x,y,on="f.eid"),updateDF)
print(merged)
