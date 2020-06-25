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
        return 0
    if a1==a and a2==a:
        return 2-2*f
    elif (a1==a and a2!=a) or (a1!=a and a2==a):
        return 1-2*f
    else:
        return -2*f
#---------------------------------------------------------------------------------------------------------------------------

# Read input files
pedfile=pd.read_table(ped,sep=" ",header=None)
pedfile=pedfile.loc[:,6:]
mapfile=pd.read_table(mapfn,names=["ID"],usecols=[1],sep="\t",header=None)
print('{:=^80}'.format(' MAP '))
print("")
print(mapfile)
print("")
mafile=pd.read_table(ma,sep="\t",header=0)
print('{:=^80}'.format(' M/A RESULTS '))
print(ma)
print(mafile)
print("")
condfile=pd.read_table(cond,header=None,names=["ID"])
print('{:=^80}'.format(' CONDITIONAL VARIANTS '))
print("")
print(condfile)
print("")

L=[]
L1=[]
# L1: variant IDs in correct order
for index, row in mapfile.iterrows():
    x=row["ID"]
    L1.append(x)
    L.append(x+"_1")
    L.append(x+"_2")

# set column names to variant id_1, id_2
pedfile.columns=L
print('{:=^80}'.format(' PED '))
print("")
print(pedfile)
print("")

# conditioning variants
L2=[]
for index, row in condfile.iterrows():
    x=row["ID"]
    L2.append(x)

var=None # variant being tested
L3=[] # conditioning variants in correct order
for x in L1:
    if x in L2:
        L3.append(x)
    else:
        var=x

print('{:=^80}'.format(' TESTED VARIANT '))
print("")
print(var)
print("")
#print(L3)

beta_var=None
f_var=None
a_var=None
tmp_betas=dict()
tmp_freqs=dict()
tmp_alleles=dict()

AL=dict() # id --> effect allele
for index, row in mafile.iterrows():
    AL[row["SNP"]]=row["A1"]
    if var==row["SNP"]:
        beta_var=float(row["b"])
        f_var=float(row["freq"])
        a_var=row["A1"]
    else:
        tmp_betas[row["SNP"]]=float(row["b"])
        tmp_freqs[row["SNP"]]=float(row["freq"])
        tmp_alleles[row["SNP"]]=row["A1"]

AF=dict() # id --> effect allele AF
for x in AL:
    t=0
    c=0
    a=AL[x]
    for i, r in pedfile.iterrows():
        a1=r[x+"_1"]
        a2=r[x+"_2"]
        if a1==a:
            c=c+1
        if a2==a:
            c=c+1
        
        if a1!="0" and a2!="0":
            t=t+1
    AF[x]=c/(2*t)

betas=[]
freqs=[]
alleles=[]

# betas etc. in correct order
for x in L2:
    betas.append(tmp_betas[x])
    freqs.append(tmp_freqs[x])
    alleles.append(tmp_alleles[x])

betas=np.asarray(betas)

print(a_var,f_var,beta_var)
print("COND ALLELES: ",alleles)
print("COND FREQS  : ",freqs)
print("COND BETAS  : ",betas)

# genotype encoded
df0=pd.DataFrame()
for index, row in mapfile.iterrows():
    x=row["ID"]
    a=mafile.loc[mafile.SNP==x,"A1"].values[0]
    f=mafile.loc[mafile.SNP==x,"freq"].values[0]
    df[x]=pedfile[[x+"_1",x+"_2"]].apply(lambda row: recode(row[0],row[1],a,AF[x]),axis=1)

print('{:=^80}'.format(' GENOTYPE MATRIX '))
print("")
print(df)
print("")

# remove rows with NaNs
#df=df0.dropna()
#print('{:=^80}'.format(' GENOTYPE MATRIX '))
#print("")
#print(df)
#print("")

#-------------------------------------------- creating necessary matrices -------------------------------------------------

# tested variant's genotypes
X2=df[[var]].to_numpy(copy=True)

# conditioning variants' genotypes
X1=df[L3].to_numpy(copy=True)

print('{:=^80}'.format(' CONDITIONING GENOTYPES (X1) '))
print("")
print(X1)
print("")
print('{:=^80}'.format(' TESTED VARIANT\'S GENOTYPES (X2) '))
print("")
print(X2)
print("")

Xp1=np.dot(np.transpose(X1),X1)
Xp2=np.dot(np.transpose(X2),X2)

print('{:=^80}'.format(' X1\'X1 '))
print("")
print(Xp1)
print("")
print('{:=^80}'.format(' X2\'X2 '))
print("")
print(Xp2)
print("")

X11=np.linalg.inv(Xp1)
X22=np.linalg.inv(Xp2)
X21=np.dot(np.transpose(X2),X1)

#print(np.dot(Xp1,X11))
#print(np.dot(Xp2,X22))

print('{:=^80}'.format(' (X1\'X1)-1 '))
print("")
print(X11)
print("")
print('{:=^80}'.format(' (X2\'X2)-1 '))
print("")
print(X22)
print("")
print('{:=^80}'.format(' (X2\'X1) '))
print("")
print(X21)
print("")

D1=np.diag(np.diagonal(Xp1))
D2=np.diag(np.diagonal(Xp2))

print('{:=^80}'.format(' D1 '))
print("")
print(D1)
print("")
print('{:=^80}'.format(' D2 '))
print("")
print(D2)
print("")

X=np.dot(np.dot(X22,np.dot(X21,X11)),D1)

print('{:=^80}'.format(' (X2\'X2)-1 X2\'X1 (X1\'X1)-1 D1 '))
print("")
print(X)
print("")

b2=beta_var-np.dot(X,betas)
print('{:=^80}'.format(' OUTPUT '))
print(ma,str(beta_var),str(b2),sep="\t")
