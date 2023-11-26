#!/usr/bin/python3.6

import argparse
import pandas as pd
import sys
import csv

def main():
    parser=argparse.ArgumentParser(description="Sort dataframe columns",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input','-i',required=True,action='store',help="Input")
    parser.add_argument('--output','-o',required=True,action='store',help="Output")
    
    if len(sys.argv[1:])==0:
        parser.print_help()
        sys.exit(0)

    try:
        args=parser.parse_args()
    except:
        sys.exit(1)

    infile=args.input
    outfile=args.output
    
    df=pd.read_table(infile,sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False)
    scols=list()
    if "f.eid" in df.columns:
        scols.append("f.eid")
    for c in sorted(df.columns):
        if c!="f.eid" and c!="CREATED" and c!="RELEASE":
            scols.append(c)
    if "CREATED" in df.columns:
        scols.append("CREATED")
    if "RELEASE" in df.columns:
        scols.append("RELEASE")
    df.to_csv(outfile,sep="\t",columns=scols,index=False,quotechar='"',quoting=csv.QUOTE_NONE,header=True)

if __name__=="__main__":
    main()
