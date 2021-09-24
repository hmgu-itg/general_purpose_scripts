#!/usr/bin/python3

import argparse
import pandas as pd
import csv
import sys
import tarfile
import os

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

C=dict()
print("Config: %s" % config)
if os.path.exists(config):
    with open(config,"r") as f:
        for line in f:
            l=line.rstrip()
            if l and not l.startswith("#"):
                x=l.split("\t",1)
                C[x[0]]=x[1]
    print(C)
    sys.exit(0)
else:
    print("ERROR: config file %s does not exist" % config,file=sys.stderr)
    sys.exit(1)
                
#-----------------------------------------------------------------------------------------------------------------------------

df=pd.read_table(infile,sep="\t",nrows=1)
L=df.columns.values.tolist()
with open(out_prefix+".txt","w") as f:
    print(*L,sep="\n",file=f)
