#!/usr/bin/python3

import sys
import os
import argparse
import re
from bisect import bisect_left

parser = argparse.ArgumentParser(description="Map GO1 variant IDs to UKBB reference panel variant IDs")
parser.add_argument('--input','-i',required=True,action="store",help="GO1 ID list")
parser.add_argument('--ukbb','-u',required=True,action="store",help="UKBB ID list")
args=parser.parse_args()

go1fname=args.input
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
    m=re.match("(\d+:\d+)_([A-Z])_([A-Z])",var)
    if m:
        return m.group(1)+"_"+m.group(3)+"_"+m.group(2)
    else:
        print("ERROR: %s is malformed" %(var),file=sys.stderr)
        return None

def findGO1InUKBBList(var,L):
    out=list()
    if binarySearch(L,var):
        out.append(var)
    var1=switchAlleles(var)
    if binarySearch(L,var1):
        out.append(var1)
    return out
#---------------------------------------------------------------------------------------------------------------------------

go1_list=list()
ukbb_list=list()

with open(go1fname) as f:
    for line in f: 
        line=line.strip()
        go1_list.append(line)

with open(ukbbfname) as f:
    ukbb_list=f.read().splitlines()

ukbb_list=sorted(ukbb_list)
    
for v in go1_list:
    L=findGO1InUKBBList(v,ukbb_list)
    if len(L)==0:
        print("# NOT FOUND %s" %(v),file=sys.stdout)
    elif len(L)==2:
        print("# AMBIGUOUS %s %s" %(v,L[0]),file=sys.stdout)
        print("# AMBIGUOUS %s %s" %(v,L[1]),file=sys.stdout)
    else:
        print("%s %s" %(v,L[0]),file=sys.stdout)
