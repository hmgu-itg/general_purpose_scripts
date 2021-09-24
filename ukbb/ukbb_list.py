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
parser.add_argument('-i','--input',required=True,action='store',help="Input file")
parser.add_argument('-o','--output',required=True,action='store',help="Output prefix")
parser.add_argument('-c','--config',required=False,action='store',help="Config file")

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    sys.exit(1)

infile=args.input
out_prefix=args.output
logF=open(out_prefix+".log","w")

if args.config is None:
    config=os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),"config.txt")
else:
    config=args.config

print("Config: %s" % config)
C=utils.readConfig(config)
if C is None:
    sys.exit(1)

#print(C)
D=utils.readDataDictionary(C["DATA_DICT"])
#print(D)

#-----------------------------------------------------------------------------------------------------------------------------

df=pd.read_table(infile,sep="\t",nrows=1)
L=df.columns.values.tolist()
H=dict()
for x in L:
    if x=="f.eid" or x=="RELEASE" or x=="CREATED":
        continue
    m=re.match("^f\.(\d+)\.\d+\.\d+",x)
    if m:
        key=m.group(1)
        H.setdefault(key,[]).append(x)
    else:
        print("ERROR: column name format error: %s" %(x),file=sys.stderr)
    
with open(out_prefix+".txt","w") as f:
    print("{}\t{}\t{}\t{}".format("Field","Instances","Description","Type"),file=f)
    for x in H:
        if x in D:
            print("{}\t{}\t{}\t{}".format(x,len(H[x]),D[x]["Field"],D[x]["ValueType"]),file=f)
        else:
            print("{}\t{}\t{}\t{}".format(x,len(H[x]),"NA","NA"),file=f)
            print("WARN: %s is not in data dictionary" % x,file=sys.stderr)
