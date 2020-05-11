#!/usr/bin/python3

import requests, sys, time
import os
import argparse
import re
import time
from requests.exceptions import Timeout
import datetime

headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

def getResponse2(server,ext,rs,headers,timeout=None):
    try:
        r = requests.get(server+ext+rs+"?",headers=headers,timeout=timeout)
        while not r.ok:
            print(str(datetime.datetime.now())+" : "+str(r.status_code)+" occured. Trying again",file=sys.stderr)
            sys.stderr.flush()
            time.sleep(10)
            r = requests.get(server+ext+rs+"?",headers=headers,timeout=timeout)
        return r.json()
    except Timeout as ex:
        print(str(datetime.datetime.now())+" : Timeout event occured", file=sys.stderr)
        sys.stderr.flush()
        return None


def parseSPDI(string):
    L=string.rsplit(":")
    #print(len(L))
    #print(len(L[3]))
    c="NA"
    m=re.search("NC_0+([^\.0]+)\.\d+",L[0])
    if m:
        c=m.group(1)
    pos=int(L[1])
    ref=L[2]
    alt=L[3]
    if len(ref)==1 and len(alt)==1:
        pos=pos+1
    return {"chr":m.group(1),"pos":pos,"ref":ref,"alt":alt}

#----------------------------------------------------------------------------------------------------------------------------------

build="grch38"

parser = argparse.ArgumentParser(description="Get chromosome, position and alleles for given rsID")
parser.add_argument('--build','-b', action="store",help="Genome build: default: grch38", default="grch38")
parser.add_argument('--rs','-r', action="store",help="rsID")
args=parser.parse_args()

if args.build!=None:
    build=args.build

rsID=args.rs
    
ext = "/variant_recoder/homo_sapiens/"
server = "http://"+build+".rest.ensembl.org"
if build=="grch38":
    server = "http://rest.ensembl.org"


#print("Current build: "+build,file=sys.stderr)

#---------------------------------------------------------------------------------------------------------------------------


r=getResponse2(server,ext,rsID,headers)
H={}
if len(r)>1:
    print("More than 1 hash for "+rsID,file=sys.stderr,flush=True)
else:
    x=r[0]
    r1=x["id"][0]
    inp=x["input"]
    if r1!=inp:
        print("WARNING: INPUT=",x["input"],"ID=",r1,file=sys.stderr,flush=True)
    else:
        H[r1]=[]
        spdi=x["spdi"]
        print(spdi)
        for z in spdi:
            m=re.search("^NC_0+",z)
            if m:
                p=parseSPDI(z)
                H[r1].append(p)

for k,r in H.items():
    print(k,",".join(list(sorted(set(x["chr"] for x in r)))),",".join(list(sorted(set(str(x["pos"]) for x in r)))),",".join(list(sorted(set(x["ref"] for x in r)))),",".join(list(sorted(set(x["alt"] for x in r)))),sep="\t")
