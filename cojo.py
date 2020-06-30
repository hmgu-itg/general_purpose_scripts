#!/usr/bin/python3

import requests, sys, time
import os
import argparse
import re
import time
import datetime
import pandas as pd
import numpy as np
from math import sqrt

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
tmp_SEs=dict()

AL=dict() # id --> effect allele
allSE=dict() # id --> SE
allbeta=dict() # id --> peta
for index, row in mafile.iterrows():
    AL[row["SNP"]]=row["A1"]
    allSE[row["SNP"]]=row["se"]
    allbeta[row["SNP"]]=row["b"]

    if var==row["SNP"]:
        beta_var=float(row["b"])
        f_var=float(row["freq"])
        a_var=row["A1"]
    else:
        tmp_betas[row["SNP"]]=float(row["b"])
        tmp_freqs[row["SNP"]]=float(row["freq"])
        tmp_alleles[row["SNP"]]=row["A1"]
        tmp_SEs[row["SNP"]]=row["se"]

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

# data of the conditioning variants
SEs=[]
betas=[]
freqs=[]
alleles=[]

# betas etc. in correct order
for x in L3:
    SEs.append(tmp_SEs[x])
    betas.append(tmp_betas[x])
    freqs.append(tmp_freqs[x])
    alleles.append(tmp_alleles[x])

betas=np.asarray(betas)

print(a_var,f_var,beta_var)
print("COND ALLELES: ",alleles)
print("COND FREQS  : ",freqs)
print("COND BETAS  : ",betas)
print("COND SE     : ",SEs)

# genotype encoded
#df0=pd.DataFrame()
df=pd.DataFrame()
for index, row in mapfile.iterrows():
    x=row["ID"]
    a=mafile.loc[mafile.SNP==x,"A1"].values[0]
    f=mafile.loc[mafile.SNP==x,"freq"].values[0]
    df[x]=pedfile[[x+"_1",x+"_2"]].apply(lambda row: recode(row[0],row[1],a,AF[x]),axis=1)

print('{:=^80}'.format(' GENOTYPE MATRIX '))
print("")
print(df)
print("")

nsamples=df.shape[0]
print('{:=^80}'.format(' N SAMPLES '))
print("")
print(nsamples)
print("")

# remove rows with NaNs
#df=df.dropna()
#print('{:=^80}'.format(' GENOTYPE MATRIX '))
#print("")
#print(df)
#print("")

#-------------------------------------------- creating necessary matrices -------------------------------------------------

allD=dict()
for v in L1:
    tmp=df[[v]].to_numpy(copy=True)
    allD[v]=np.dot(np.transpose(tmp),tmp)
    
print('{:=^80}'.format(' D '))
print("")
for v in L1:
    print(v,allD[v],sep="\t")
print("")

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

#------------------------------------------------------ calculating conditional beta ----------------------------------------------

b2=beta_var-np.dot(X,betas)
print('{:=^80}'.format(' OUTPUT '))
print(ma,str(beta_var),str(b2),sep="\t")

#-------------------------------------------------------------calculating SE ------------------------------------------------------
tmpL=[]
for v in L1:
    tmpL.append((nsamples-1)*allSE[v]*allSE[v]+allD[v]*allbeta[v]*allbeta[v])

yty=np.median(tmpL)
#yty=np.mean(tmpL)
print('{:=^80}'.format(' yty '))
print("")
print(tmpL)
print(yty)
print("")

tmp=np.dot(X11,D1)
b1=np.dot(tmp,betas)
print('{:=^80}'.format(' b1 '))
print("")
print(b1)
print("")

tmp=np.dot(np.transpose(b1),D1)
a2=np.dot(tmp,betas)
print('{:=^80}'.format(' a2 '))
print("")
print(a2)
print("")

sigma2=(yty-a2-b2*allD[var]*allbeta[var])/(nsamples-len(L1))
print('{:=^80}'.format(' sigma2 '))
print("")
print(sigma2)
print("")

t1=np.dot(np.transpose(X2),X1)
t2=np.dot(np.transpose(X1),X2)
z=D2-np.dot(np.dot(t1,X11),t2)
print('{:=^80}'.format(' z '))
print("")
print(z)
print("")

se2=sqrt(sigma2*z)/allD[var]
print('{:=^80}'.format(' SE2 '))
print("")
print(se2)
print("")
