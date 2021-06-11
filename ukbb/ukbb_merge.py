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
parser.add_argument("-k", "--keep",required=False,action='store_true',help="Also include new IDs")
if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    sys.exit(1)

mode="right"
if args.keep:
    mode="outer"

infile=args.input
updates=args.update
out_prefix=args.output
logF=open(out_prefix+".log","w")
datestr=datetime.datetime.now().strftime("%F")

print("input: %s\nupdates: %s\noutput prefix: %s\n" %(infile,",".join(updates),out_prefix),file=logF)

inputDF=pd.read_table(infile,sep="\t",header=[0,1],dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False)
print("Input rows: %d" % len(inputDF),file=logF)
print("Input columns: %d" % len(inputDF.columns.values.tolist()),file=logF)

input_classes=dict()
for x in inputDF.columns.values.tolist():
    if not x[0] in ("CREATED","RELEASE","f.eid"):
        input_classes[x[0]]=x[1]
print("Input column classes: %d" % len(set(input_classes.values())),file=logF)
for v in sorted(list(set(input_classes.values()))):
    L=[x for x in input_classes.keys() if input_classes[x]==v]
    print("Column class %s: %d column(s)" % (v,len(L)),file=logF)
print("",file=logF)        
    
release=int(inputDF["RELEASE"].iloc[2])
print("previous release: %d\n" % release,file=logF)
release+=1
        
updateDF=list()
for f in updates:
    updateDF.append(pd.read_table(f,sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False))

print("Checking IDs in update DFs\n",file=logF)
for i in range(len(updateDF)):
    print("Update %d rows: %d" % (i,len(updateDF[i])),file=logF)
    print("Update %d columns: %d\n" % (i,len(updateDF[i].columns.values.tolist())),file=logF)
    for j in range(len(updateDF)):
        if i<j:
            if not updateDF[i]["f.eid"].equals(updateDF[j]["f.eid"]):
                print("ERROR: IDs in %d,%d not equal" %(i,j),file=logF)
                sys.exit(1)
print("Done\n",file=logF)

print("Checking if columns in update DFs are disjoint",file=logF)
for i in range(len(updateDF)):
    for j in range(len(updateDF)):
        if i<j:
            L=[x for x in updateDF[i].columns.values.tolist() if x in updateDF[j].columns.values.tolist()]
            if len(L)>1:
                print("ERROR: %d,%d have common columns:" %(i,j),file=logF)
                for c in L:
                    if c!="f.eid":
                        print("%s" % c,file=logF)
                print("",file=logF)        
                sys.exit(1)
print("Done\n",file=logF)

update_classes=dict()
cur_c=0
for df in updateDF:
    for c in df.columns.values.tolist():
        if c=="f.eid":
            continue
        update_classes[c]=cur_c
    cur_c+=1
maxc=max(update_classes.values())

print("Merging update DFs",file=logF)
if len(updateDF)==1:
    merged=updateDF[0]
else:
    merged=reduce(lambda x, y:pd.merge(x,y,on="f.eid"),updateDF)
print("Done\n",file=logF)
affected_classes=set()
for c in merged.columns.values.tolist():
    if c in input_classes:
        affected_classes.add(input_classes[c])
columns_to_remove=list(["CREATED","RELEASE"])
for c in input_classes:
    if c=="f.eid":
        continue
    if input_classes[c] in affected_classes:
        columns_to_remove.append(c)
inputDF.columns=[x[0] for x in inputDF.columns.values.tolist()]
print("Common IDs between input and update: %d" %(len(set(inputDF["f.eid"]).intersection(set(merged["f.eid"])))),file=logF)
print("IDs in input but not in update (will be removed): %d" %(len(set(inputDF["f.eid"]).difference(set(merged["f.eid"])))),file=logF)
print("IDs in update but not in input: %d" %(len(set(merged["f.eid"]).difference(set(inputDF["f.eid"])))),file=logF)
df=inputDF.drop(columns=columns_to_remove)
# keep new IDs in, remove IDs not present in the update
df2=pd.merge(df,merged,on="f.eid",how=mode)
df2.replace(np.nan,"NA",inplace=True)
df2["RELEASE"]=release
df2["CREATED"]=datestr
new_classes=dict()
temp=dict()
for c in df2.columns.values.tolist():
    if c in ["f.eid","RELEASE","CREATED"]:
        continue
    else:
        if c in update_classes:
            new_classes[c]=update_classes[c]
        else:
            c1=input_classes[c]
            if c1 not in temp:
                temp[c1]=maxc+1
                maxc+=1
            new_classes[c]=temp[c1]
print("Output column classes: %d" % len(set(new_classes.values())),file=logF)
for v in sorted(list(set(new_classes.values()))):
    L=[x for x in new_classes.keys() if new_classes[x]==v]
    print("Column class %s: %d column(s)" % (v,len(L)),file=logF)
print("",file=logF)        
for c in ["f.eid","RELEASE","CREATED"]:
    new_classes[c]="NA"
L=[(x,new_classes[x]) for x in df2.columns.values.tolist()]
df2.columns=pd.MultiIndex.from_tuples(L)
print("Output rows: %d" % len(df2),file=logF)
print("Output columns: %d" % len(df2.columns.values.tolist()),file=logF)
df2.to_csv(out_prefix+".txt.gz",sep="\t",index=False)
