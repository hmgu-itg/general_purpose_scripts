#!/usr/bin/python3

import argparse
import pandas as pd
import numpy as np
import csv
import sys
import datetime
from functools import reduce

parser=argparse.ArgumentParser()
parser.add_argument('-i','--input',required=True,action='store',help="Input file")
parser.add_argument('-u','--update',required=True,action='append',help="Update file(s)")
parser.add_argument('-o','--output',required=True,action='store',help="Output prefix")

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    sys.exit(1)

infile=args.input
updates=args.update
out_prefix=args.output

print("input: %s\nupdates: %s\n" %(infile,",".join(updates)))

inputDF=pd.read_table(infile,sep="\t",header=[0,1],dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False)
#print(inputDF)
#print(inputDF.columns.values.tolist())
input_classes=dict()
for x in inputDF.columns.values.tolist():
    if not x[0] in ("CREATED","RELEASE","f.eid"):
        print("%s\t%s" % (x[0],x[1]))
        input_classes[x[0]]=x[1]

release=int(inputDF["RELEASE"].iloc[1])
print("release=%d" % release)
release+=1
        
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


update_classes=dict()
cur_c=0
for df in updateDF:
    for c in df.columns.values.tolist():
        if c=="f.eid":
            continue
        update_classes[c]=cur_c
    cur_c+=1
maxc=max(update_classes.values())

print(update_classes)
print("Merging update DFs")
merged=reduce(lambda x, y:pd.merge(x,y,on="f.eid"),updateDF)
print(merged)
affected_classes=set()
for c in merged.columns.values.tolist():
    if c in input_classes:
        affected_classes.add(input_classes[c])
print(affected_classes)
columns_to_remove=list(["CREATED","RELEASE"])
for c in input_classes:
    if c=="f.eid":
        continue
    if input_classes[c] in affected_classes:
        columns_to_remove.append(c)
print(columns_to_remove)
#print(inputDF)
inputDF.columns=[x[0] for x in inputDF.columns.values.tolist()]
df=inputDF.drop(columns=columns_to_remove)
#print(inputDF)
datestr=datetime.datetime.now().strftime("%F")
print(df)
df2=pd.merge(df,merged,on="f.eid",how="outer")
df2.replace(np.nan,"NA",inplace=True)
df2["RELEASE"]=release
df2["CREATED"]=datestr
print(df2)
new_classes=dict()
temp=dict()
for c in df2.columns.values.tolist():
    if c in ["f.eid","RELEASE","CREATED"]:
        new_classes[c]="NA"
    else:
        if c in update_classes:
            new_classes[c]=update_classes[c]
        else:
            c1=input_classes[c]
            if c1 not in temp:
                temp[c1]=maxc+1
                maxc+=1
            new_classes[c]=temp[c1]
L=[(x,new_classes[x]) for x in df2.columns.values.tolist()]
df2.columns=pd.MultiIndex.from_tuples(L)
print(df2)
df2.to_csv(out_prefix+".txt",sep="\t",index=False)
