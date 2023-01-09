#!/usr/bin/python3

import sys
import os
import argparse
import re
from bisect import bisect_left

parser = argparse.ArgumentParser(description="Map GO2 variant IDs to 1KG reference panel variant IDs")
parser.add_argument('--input','-i',required=True,action="store",help="GO2 ID list")
parser.add_argument('--g','-g',required=True,action="store",help="1KG ID list")
args=parser.parse_args()

go2fname=args.input
kfname=args.g

#---------------------------------------------------------------------------------------------------------------------------
# if sorted L contains x
def binarySearch(L,x): 
    i=bisect_left(L,x) 
    if i!=len(L) and L[i]==x: 
        return True
    else: 
        return False

def switchAlleles(var):
    m=re.match("(\d+:\d+)_([A-Z]+)_([A-Z]+)$",var)
    if m:
        return m.group(1)+"_"+m.group(3)+"_"+m.group(2)
    else:
        print("ERROR: %s is malformed" %(var),file=sys.stderr)
        return None

def findGO2In1KGList(var,L):
    out=list()
    if isSNP(var):
        if binarySearch(L,var):
            out.append(var)
        var1=switchAlleles(var)
        if binarySearch(L,var1):
            out.append(var1)
    else:
        var1,var2=recodeIndel(var)
        # print("%s,%s %s" %(var1,var2,var),file=sys.stderr)
        if binarySearch(L,var1):
            out.append(var1)
        if binarySearch(L,var2):
            out.append(var2)
    return out

# for 1:12345_A_ACG retruns 1:12345_D_I,1:12345_I_D
def recodeIndel(var):
    m=re.match("(\d+:\d+)_([ACGT]+)_([ACGT]+)",var)
    if m:
        return m.group(1)+"_D_I",m.group(1)+"_I_D"
    else:
        print("ERROR: %s is malformed" %(var),file=sys.stderr)
        return None
    
def isSNP(var):
    return re.match("\d+:\d+_[ACGT]_[ACGT]$",var)

#---------------------------------------------------------------------------------------------------------------------------

go2_list=list()
k_list=list()

with open(go2fname) as f:
    for line in f: 
        line=line.strip()
        go2_list.append(line)

k_new2old=dict()
with open(kfname) as f:
    for line in f: 
        line=line.strip()
        L=line.split()
        L1=L[1].split(sep=":")
        if len(L1)==4:
            if re.match("rs\d+",L1[0]):
                new_id=L[0]+":"+L1[1]+"_"+L1[2]+"_"+L1[3]
            else:
                new_id=L1[0]+":"+L1[1]+"_"+L1[2]+"_"+L1[3]
            k_list.append(new_id)
            k_new2old[new_id]=L[1]
        else:
            print("SKIPPING %s" %(L[1]),file=sys.stderr)

k_list=sorted(k_list)

for v in go2_list:
    if binarySearch(k_list,v):
        print("%s %s" %(v,k_new2old[v]),file=sys.stdout)
    else:
        print("# NOT FOUND %s" %(v),file=sys.stdout)
