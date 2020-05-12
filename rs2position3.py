#!/usr/bin/python3

import requests, sys, time
import os
import argparse
import re
import time
from requests.exceptions import Timeout
import datetime

headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

def list2string(snps):
    return "{\"ids\":["+",".join(snps)+"]}"

def getResponse2(request_string,headers,data,timeout=None,max_attempts=-1):
    attempt=1
    try:
        r = requests.post(request_string,headers=headers,data=data,timeout=timeout)
        #print(data)
        while not r.ok:
            time.sleep(10)
            print(str(datetime.datetime.now())+" : Error "+str(r.status_code)+" occured. Trying again",file=sys.stderr)
            sys.stderr.flush()
            attempt+=1
            if max_attempts!=-1 and attempt==max_attempts:
                return None

            r = requests.post(request_string,headers=headers,data=data)
        return r.json()
    except Timeout as ex:
        print(str(datetime.datetime.now())+" : Timeout event occured", file=sys.stderr)
        sys.stderr.flush()
        return None

def parseSPDI(string):
    L=string.rsplit(":")
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

def printHash(H):
    for k,r in H.items():
        print(k,",".join(list(sorted(set(x["chr"] for x in r)))),",".join(list(sorted(set(str(x["pos"]) for x in r)))),",".join(list(sorted(set(x["ref"] for x in r)))),",".join(list(sorted(set(x["alt"] for x in r)))),sep="\t",flush=True)

#----------------------------------------------------------------------------------------------------------------------------------

build="grch38"
batchsize=200

parser = argparse.ArgumentParser(description="Get chromosome, position, REF and ALT alleles for a list of rsIDs\nINPUT: STDIN\nOUTPUT: STDOUT\n./rs2position3.py <INFILE >OUTFILE")
parser.add_argument('--build','-b', action="store",help="Genome build: default: grch38", default="grch38")
parser.add_argument('--size','-s', action="store",help="Batch size: default: 200", default=200)
args=parser.parse_args()

if args.build!=None:
    build=args.build

if args.size!=None:
    batchsize=int(args.size)

ext = "/variant_recoder/homo_sapiens"
server = "https://"+build+".rest.ensembl.org"
if build=="grch38":
    server = "https://rest.ensembl.org"

#print("Current build: "+build,file=sys.stderr)

#---------------------------------------------------------------------------------------------------------------------------

IDs=sys.stdin.readlines()
#print("Read "+str(len(IDs))+" IDs",file=sys.stderr)

total_lines=len(IDs)

cur_line=0
while cur_line<total_lines:
    L=[]
    for i in range(cur_line,min(cur_line+batchsize,total_lines)):
        L.append("\""+IDs[i].rstrip(os.linesep)+"\"")    
    r=getResponse2(server+ext,headers,list2string(L),max_attempts=5)
    if r:
        for snprec in r:
            H={}
            id1=snprec["id"][0]
            H[id1]=[]
            spdi=snprec["spdi"]
            for z in spdi:
                m=re.search("^NC_0+",z)
                if m:
                    p=parseSPDI(z)
                    H[id1].append(p)
            printHash(H)
    else:
        for i in range(cur_line,min(cur_line+batchsize,total_lines)):
            print(IDs[i].rstrip(os.linesep),"NA","NA","NA",file=sys.stdout,flush=True)    
                    
    cur_line+=batchsize
