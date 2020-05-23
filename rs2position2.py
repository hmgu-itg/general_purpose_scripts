#!/usr/bin/python3

import requests, sys, time
import os
import argparse
import re
import time
from requests.exceptions import Timeout
import datetime

headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

def getResponse2(server,ext,rs,headers,timeout=None,max_attempts=-1):
    attempt=1
    try:
        r = requests.get(server+ext+rs+"?",headers=headers,timeout=timeout)
        while not r.ok:
            time.sleep(10)
            print(str(datetime.datetime.now())+" : Error "+str(r.status_code)+" occured. Trying again",file=sys.stderr)
            sys.stderr.flush()
            attempt+=1
            if max_attempts!=-1 and attempt==max_attempts:
                return None

            r = requests.get(server+ext+rs+"?",headers=headers,timeout=timeout)
            try:
                ret=r.json()
                return ret
            except ValueError:
                print(str(datetime.datetime.now())+" : JSON decoding error", file=sys.stderr)
                sys.stderr.flush()
                return None
    except Timeout as ex:
        print(str(datetime.datetime.now())+" : Timeout exception occured", file=sys.stderr)
        sys.stderr.flush()
        return None
    except TooManyRedirects as ex:
        print(str(datetime.datetime.now())+" : TooManyRedirects exception occured", file=sys.stderr)
        sys.stderr.flush()
        return None
    except RequestException as ex:
        print(str(datetime.datetime.now())+" : RequestException occured", file=sys.stderr)
        sys.stderr.flush()
        return None

def parseSPDI(string):
    L=string.rsplit(":")
    #print(len(L))
    #print(len(L[3]))
    c=L[0]
    m=re.search("NC_0+(\d+)\.\d+",L[0])
    if m:
        c=m.group(1)
    pos=int(L[1])
    ref=L[2]
    alt=L[3]
    if len(ref)==1 and len(alt)==1:
        pos=pos+1
    return {"chr":c,"pos":pos,"ref":ref,"alt":alt}

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

timeout=60
max_attempts=5

r=getResponse2(server,ext,rsID,headers,timeout,max_attempts)
H={}

if r:
    print(repr(r))
    if len(r)>1:
        print("WARNING: More than 1 hash for "+rsID,file=sys.stderr,flush=True)
        
    x=r[0]
    r1=x["id"][0]
    if r1!=rsID:
        print("WARNING: INPUT ID="+rsID,"RETRIEVED ID="+r1,file=sys.stderr,flush=True)

    H[rsID]=[]
    spdi=x["spdi"]

    for z in spdi:
        m=re.search("^NC_0+",z)
        if m:
            p=parseSPDI(z)
            H[rsID].append(p)

    s=H[rsID]
    positions=set(x["chr"]+":"+str(x["pos"]) for x in s)
    if len(positions)>1:
        print("ERROR: more than one position for "+rsID,file=sys.stderr,flush=True)
    elif len(positions)<1:
        print("ERROR: no position for "+rsID,file=sys.stderr,flush=True)
    else:
        L=positions.pop().rsplit(":")
        print(rsID,L[0],L[1],sep='\t',file=sys.stdout,flush=True)
else:
    print("ERROR: getResponse2 returned None for "+rsID,file=sys.stderr,flush=True)    
    print(rsID,"NA","NA",sep='\t')


#    print(k,",".join(list(sorted(set(x["chr"] for x in r)))),",".join(list(sorted(set(str(x["pos"]) for x in r)))),",".join(list(sorted(set(x["ref"] for x in r)))),",".join(list(sorted(set(x["alt"] for x in r)))),sep="\t")
