#!/usr/bin/python3

import requests, sys, time
import os
import argparse
import re
import time
import datetime
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description="Perform conditional analysis")
parser.add_argument('--ma','-m', action="store",help="Meta-analysis results")
parser.add_argument('--ped','-p', action="store",help="PED file")
parser.add_argument('--map','-map', action="store",help="MAP file")
parser.add_argument('--cond','-c', action="store",help="Conditioning signals")
args=parser.parse_args()

ma=args.ma
ped=args.ped
mapfn=args.map
cond=args.cond

#---------------------------------------------------------------------------------------------------------------------------
def recode(a1,a2,a,f):
    if a1=="0" or a2=="0":
        return None
    if a1==a and a2==a:
        return -2*f
    elif (a1==a and a2!=a) or (a1!=a and a2==a):
        return 1-2*f
    else:
        return 2-2*f

#---------------------------------------------------------------------------------------------------------------------------

# Read input files

pedfile=pd.read_table(ped,sep=" ",header=None)
pedfile=pedfile.loc[:,6:]
mapfile=pd.read_table(mapfn,names=["ID"],usecols=[1],sep="\t",header=None)
print(mapfile)
mafile=pd.read_table(ma,sep="\t",header=0)
print(mafile)
condfile=pd.read_table(cond,header=None,names=["ID"])
print(condfile)
print("")

S1=set(mapfile["ID"])
S2=set(condfile["ID"])

# variant being tested
var=S1.difference(S2).pop()
print(var)

beta_var=None
f_var=None
a_var=None
betas=[]
freqs=[]
alleles=[]
for index, row in mafile.iterrows():
    if var==row["SNP"]:
        beta_var=row["b"]
        f_var=row["freq"]
        a_var=row["A1"]
    else:
        betas.append(row["b"])
        freqs.append(row["freq"])
        alleles.append(row["A1"])

print(a_var,f_var,beta_var)
print(alleles,freqs,betas)

L=[]
for index, row in mapfile.iterrows():
    x=row["ID"]
    if x==var:
        beta_var=row
    L.append(x+"_1")
    L.append(x+"_2")
    f=mafile.loc[mafile["SNP"]==x,"freq"].values[0]
    a=mafile.loc[mafile["SNP"]==x,"A1"].values[0]
    print(index,row["ID"],f,a)

pedfile.columns=L
print(pedfile)
df=pd.DataFrame()
for index, row in mapfile.iterrows():
    x=row["ID"]
    f=mafile.loc[mafile["SNP"]==row["ID"],"freq"].values[0]
    a=mafile.loc[mafile["SNP"]==row["ID"],"A1"].values[0]
    df[x]=pedfile[[x+"_1",x+"_2"]].apply(lambda row: recode(row[0],row[1],a,f),axis=1)

print(df)

