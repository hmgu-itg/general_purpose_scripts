#!/usr/bin/python3

import argparse
import pandas as pd
import csv
import sys
import tarfile
import os
from itgukbb import utils
import re

parser=argparse.ArgumentParser()
parser.add_argument('--input','-i',required=True,action='store',help="Input file")
parser.add_argument('--output','-o',required=True,action='store',help="Output prefix")
parser.add_argument('--config','-c',required=False,action='store',help="Config file")
parser.add_argument('--list','-l',required=False,action='store_true',help="Output column information")
parser.add_argument('--majority','-majority',required=False,action='append',help="Output most frequent value, over all instances")
parser.add_argument('--mean','-mean',required=False,action='append',help="Output mean, over all instances")
parser.add_argument('--min-missing','-min-missing',required=False,action='append',help="Output instance with least NAs",dest="min_missing")
parser.add_argument('--all','-a',required=False,action='append',help="Output all instances")
parser.add_argument('--cc','-cc',required=False,action='append',help="Recode as case/control (1/0)")

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    sys.exit(1)

to_list=args.list
infile=args.input
out_prefix=args.output
logF=open(out_prefix+".log","w")

#-----------------------------------------------------------------------------------------------------------------------------

if args.config is None:
    config=os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),"config.txt")
else:
    config=args.config

print("Config: %s" % config)
C=utils.readConfig(config)
if C is None:
    sys.exit(1)

# key: field
D=utils.readDataDictionary(C["DATA_DICT"])

df=pd.read_table(infile,sep="\t",nrows=1)
L=df.columns.values.tolist()
H=dict() # short field name --> list of matching full names
for x in L:
    if x=="f.eid" or x=="RELEASE" or x=="CREATED":
        continue
    m=re.match("^f\.(\d+)\.\d+\.\d+$",x)
    if m:
        key=m.group(1)
        H.setdefault(key,[]).append(x)
    else:
        print("ERROR: column name format error: %s" %(x),file=sys.stderr)

#-----------------------------------------------------------------------------------------------------------------------------

if to_list:
    with open(out_prefix+".txt","w") as f:
        print("{}\t{}\t{}\t{}".format("Field","Instances","Description","Type"),file=f)
        for x in H:
            if x in D:
                print("{}\t{}\t{}\t{}".format(x,len(H[x]),D[x]["Field"],D[x]["ValueType"]),file=f)
            else:
                print("{}\t{}\t{}\t{}".format(x,len(H[x]),"NA","NA"),file=f)
                print("WARN: %s is not in data dictionary" % x,file=sys.stderr)
else:
    print("INFO: selecting fields",file=sys.stderr)
    ccfields=dict() # field --> set of "case" values
    if not args.cc is None:
        for x in args.cc:
            a=x.split(":",1)
            if len(a)!=2:
                print("ERROR: wrong --cc option value: %s" % x)
                continue
            b=a[1].split(",")
            for z in b:
                ccfields.setdefault(a[0],set()).add(z)
    majfields=set()
    if not args.majority is None:
        for x in args.majority:
            majfields.add(x)
    meanfields=set()
    if not args.mean is None:
        for x in args.mean:
            meanfields.add(x)
    minmissfields=set()
    if not args.min_missing is None:
        for x in args.min_missing:
            minmissfields.add(x)
    allfields=set()
    if not args.all is None:
        for x in args.all:
            allfields.add(x)

    # check field availability
    tmp=set()
    for s in [ccfields,majfields,meanfields,minmissfields,allfields]:
        for f in s:
            if not f in H:
                print("WARN: field %s is not in input header" % f,file=sys.stderr)
                tmp.add(f)
        for f in tmp:
            if s==ccfields:
                del s[f]
            else:
                s.remove(f)
        tmp.clear()

    # check field type
    tmp.clear()
    for f in ccfields:
        t=D[f]["ValueType"]
        if t!="Categorical single" and t!="Categorical multiple":
            print("WARN: case/control field %s has type %s; only \"Categorical single\" or \"Categorical multiple\" type is allowed for case/control" % (f,t),file=sys.stderr)
            tmp.add(f)
    for x in tmp:
        del ccfields[x]
    tmp.clear()

    for f in majfields:
        t=D[f]["ValueType"]
        if t!="Categorical single":
            print("WARN: majority field %s has type %s; only \"Categorical single\" type is allowed for majority" % (f,t),file=sys.stderr)
            tmp.add(f)
    for x in tmp:
        del majfields[x]
    tmp.clear()

    for f in meanfields:
        t=D[f]["ValueType"]
        if t!="Continuous" and t!="Integer":
            print("WARN: mean field %s has type %s; only \"Continuous\" and \"Integer\" types are allowed for mean" % (f,t),file=sys.stderr)
            tmp.add(f)
    for x in tmp:
        del meanfields[x]
    tmp.clear()

    to_keep=set()
    to_keep.add("f.eid")
    for s in [ccfields,majfields,meanfields,minmissfields,allfields]:
        for f in s:
            for x in H[f]:
                to_keep.add(x)

    print("INFO: reading columns %s" % repr(to_keep),file=sys.stderr)
    df=pd.read_table(infile,skiprows=[1],sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=to_keep)
    df.to_csv(out_prefix+".txt",sep="\t",index=False,quotechar='"',quoting=csv.QUOTE_NONE)
