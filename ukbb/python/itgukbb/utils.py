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

def customMean(L):
    s=0.0
    c=0
    for x in L:
        if x!="NA":
            s+=x
            c+=1
    if c==0:
        return "NA"
    else:
        return s/c

def customMajority(L):
    H=dict()
    f=1
    for x in L:
        if x!="NA":
            f=0
            H[x]=H.setdefault(x,0)+1
    if f==1:
        return "NA"
    else:
        return max(H,key=H.get)

def customCC(L,values):
    f=1
    for x in L:
        if x!="NA":
            f=0
            if x in values:
                return 1
    if f==1:
        return "NA"
    else:
        return 0

def addSummaryColumn(df,columns,new_name,values,method):
    if method=="mean":
        df[new_name]=df[columns].agg(customMean,axis="columns")
    elif method=="majority":
        df[new_name]=df[columns].agg(customMajority,axis="columns")
    elif method=="cc":
        df[new_name]=df[columns].agg(lambda L:customCC(L,values),axis="columns")
    elif method=="minmissing":
        H=dict()
        for c in columns:
            H[c]=(df[c]=="NA").sum()
        c0=min(H,key=H.get)
        df[new_name]=df[c0]
    else:
        print("ERROR: unknown method %s" % (method))
    return df
