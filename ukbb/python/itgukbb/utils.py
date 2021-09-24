import os
import sys
import pandas as pd
import csv

def readConfig(fname):
    C=dict()
    if os.path.exists(fname):
        with open(fname,"r") as f:
            for line in f:
                l=line.rstrip()
                if l and not l.startswith("#"):
                    x=l.split("\t",1)
                    C[x[0]]=x[1]
        return C
    else:
        print("ERROR: config file %s does not exist" % fname,file=sys.stderr)
        return None

def readDataDictionary(fname):
    if os.path.exists(fname):
        df=pd.read_table(fname,sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["FieldID","Field","ValueType"])
        return df.set_index("FieldID").T.to_dict()
    else:
        print("ERROR: dictionary file %s does not exist" % fname,file=sys.stderr)
        return None
