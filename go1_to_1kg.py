#!/usr/bin/python3

import sys
import os
import argparse
import re
from bisect import bisect_left

parser = argparse.ArgumentParser(description="Map GO1 variant IDs to 1KG reference panel variant IDs")
parser.add_argument('--input','-i',required=True,action="store",help="GO1 ID list")
parser.add_argument('--g','-g',required=True,action="store",help="1KG ID list")
args=parser.parse_args()

go1fname=args.input
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

def findGO1In1KGList(var,L):
    out=list()
    if binarySearch(L,var):
        out.append(var)
    var1=switchAlleles(var)
    if binarySearch(L,var1):
        out.append(var1)
    return out
#---------------------------------------------------------------------------------------------------------------------------

go1_list=list()
k_list=list()

with open(go1fname) as f:
    for line in f: 
        line=line.strip()
        go1_list.append(line)

k_new2old=dict()
with open(kfname) as f:
    for line in f: 
        line=line.strip()
        L=line.split()
        L1=L[1].split(sep=":")
        if len(L1)==4:
            if len(L1[2])>len(L1[3]):
                L1[2]="I"
                L1[3]="D"
            elif len(L1[2])<len(L1[3]):
                L1[2]="D"
                L1[3]="I"                
            if re.match("rs\d+",L1[0]):
                new_id=L[0]+":"+L1[1]+"_"+L1[2]+"_"+L1[3]
            else:
                new_id=L1[0]+":"+L1[1]+"_"+L1[2]+"_"+L1[3]
            k_list.append(new_id)
            k_new2old[new_id]=L[1]
        else:
            print("SKIPPING %s" %(L[1]),file=sys.stderr)

k_list=sorted(k_list)
    
for v in go1_list:
    L=findGO1In1KGList(v,k_list)
    if len(L)==0:
        print("# NOT FOUND %s" %(v),file=sys.stdout)
    elif len(L)==2:
        print("# AMBIGUOUS %s %s" %(v,L[0]),file=sys.stdout)
        print("# AMBIGUOUS %s %s" %(v,L[1]),file=sys.stdout)
    else:
        print("%s %s" %(v,k_new2old[L[0]]),file=sys.stdout)
