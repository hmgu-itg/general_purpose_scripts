#!/usr/bin/python3

import sys
import os
import argparse
import re
from bisect import bisect_left

parser = argparse.ArgumentParser(description="Map GO2 variant IDs to UKBB reference panel variant IDs")
parser.add_argument('--input','-i',required=True,action="store",help="GO2 ID list")
parser.add_argument('--ukbb','-u',required=True,action="store",help="UKBB ID list")
args=parser.parse_args()

go2fname=args.input
ukbbfname=args.ukbb

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

def findGO2InUKBBList(var,L):
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
ukbb_list=list()

with open(go2fname) as f:
    for line in f: 
        line=line.strip()
        go2_list.append(line)

with open(ukbbfname) as f:
    ukbb_list=f.read().splitlines()

go2_list=sorted(go2_list)
ukbb_list=sorted(ukbb_list)

dups=set()
for v in go2_list:
    v2=switchAlleles(v)
    if binarySearch(go2_list,v2):
        dups.add(v)
        dups.add(v2)

for v in go2_list:
    if v in dups:
        print("# DUPLICATE %s" %(v),file=sys.stdout)
    else:
        L=findGO2InUKBBList(v,ukbb_list)
        if len(L)==0:
            print("# NOT FOUND %s" %(v),file=sys.stdout)
        elif len(L)==2:
            print("# AMBIGUOUS %s %s" %(v,L[0]),file=sys.stdout)
            print("# AMBIGUOUS %s %s" %(v,L[1]),file=sys.stdout)
        else:
            print("%s %s" %(v,L[0]),file=sys.stdout)
