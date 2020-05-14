#!/usr/bin/python3

import requests, sys, time
import os
import argparse
import re
import time
from requests.exceptions import Timeout
import datetime

headers={ "Content-Type" : "application/json"}

def getResponse2(server,ext,headers,timeout=None,max_attempts=-1):
    attempt=1
    try:
        r = requests.get(server+ext,headers=headers,timeout=timeout)
        while not r.ok:
            print(str(datetime.datetime.now())+" : "+str(r.status_code)+" occured. Trying again",file=sys.stderr)
            sys.stderr.flush()
            attempt+=1
            if max_attempts!=-1 and attempt>max_attempts:
                return None
            
            time.sleep(10)
            r = requests.get(server+ext,headers=headers,timeout=timeout)
        
        return r.json()
    except Timeout as ex:
        print(str(datetime.datetime.now())+" : Timeout event occured", file=sys.stderr)
        sys.stderr.flush()
        return None

#----------------------------------------------------------------------------------------------------------------------------------

build="grch38"
a1=None
a2=None

parser = argparse.ArgumentParser(description="Get rs ID(s) for given chromosome, position and (optionally) alleles")
parser.add_argument('--build','-b', action="store",help="Genome build: default: grch38", default="grch38")
parser.add_argument('--chrom','-c', action="store",help="chromosome")
parser.add_argument('--pos','-p', action="store",help="position")
parser.add_argument('--a1','-a1', action="store",help="allele 1")
parser.add_argument('--a2','-a2', action="store",help="allele 2")
args=parser.parse_args()

if args.build!=None:
    build=args.build

chrom=args.chrom
pos=int(args.pos)
alleles0=[]
if args.a1!=None:
    a1=args.a1
    alleles0.append(a1)
else:
    alleles0.append("NA")

if args.a2!=None:
    a2=args.a2
    alleles0.append(a2)
else:
    alleles0.append("NA")

alleles0.sort()

#print(*alleles0,sep="\t");

# default: exact position (for SNPs)    
ext = "/overlap/region/human/"+chrom+":"+str(pos)+"-"+str(pos)+"?feature=variation"
server = "http://"+build+".rest.ensembl.org"
if build=="grch38":
    server = "http://rest.ensembl.org"


#print("Current build: "+build,file=sys.stderr)

#---------------------------------------------------------------------------------------------------------------------------

timeout=60
max_attempts=5

# TODO: fix
if a1 and a2 and (len(a1)!=1 or len(a2)!=1):
    ext = "/overlap/region/human/"+chrom+":"+str(pos-200)+"-"+str(pos+200)+"?feature=variation"
#

r=getResponse2(server,ext,headers,timeout,max_attempts)

H={}

if r:
    #print(repr(r))

    for x in r:
        alleles=x["alleles"]
        c=x["seq_region_name"]
        strand=x["strand"]
        rsID=x["id"]
        start=int(x["start"])
        end=int(x["end"])
        #if start==pos+1:
        alleles.sort()
        print(rsID,c,str(start),str(end),",".join(alleles),chrom,str(pos),",".join(alleles0),sep="\t")
        

#     r1=x["id"][0]
#     if r1!=rsID:
#         print("WARNING: INPUT ID="+rsID,"RETRIEVED ID="+r1,file=sys.stderr,flush=True)

#     H[rsID]=[]
#     spdi=x["spdi"]

#     for z in spdi:
#         m=re.search("^NC_0+",z)
#         if m:
#             p=parseSPDI(z)
#             H[rsID].append(p)

#     s=H[rsID]
#     positions=set(x["chr"]+":"+str(x["pos"]) for x in s)
#     if len(positions)>1:
#         print("ERROR: more than one position for "+rsID,file=sys.stderr,flush=True)
#     elif len(positions)<1:
#         print("ERROR: no position for "+rsID,file=sys.stderr,flush=True)
#     else:
#         L=positions.pop().rsplit(":")
#         print(rsID,L[0],L[1],sep='\t',file=sys.stdout,flush=True)
# else:
#     print("ERROR: getResponse2 returned None for "+rsID,file=sys.stderr,flush=True)    
#     print(rsID,"NA","NA",sep='\t')


#    print(k,",".join(list(sorted(set(x["chr"] for x in r)))),",".join(list(sorted(set(str(x["pos"]) for x in r)))),",".join(list(sorted(set(x["ref"] for x in r)))),",".join(list(sorted(set(x["alt"] for x in r)))),sep="\t")
