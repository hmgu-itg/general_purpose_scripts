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

def getResponse(request_string,headers,data,max_attempts=-1):
    attempt=1
    r = requests.post(request_string,headers=headers,data=data)
    while not r.ok:
        #print("attempt "+str(attempt),file=sys.stderr)
        attempt+=1
        if max_attempts!=-1 and attempt==max_attempts:
            return None
        time.sleep(15)
        r = requests.post(request_string,headers=headers,data=data)
    return r.json()

def getResponse2(request_string,headers,data,timeout=None):
    try:
        r = requests.post(request_string,headers=headers,data=data,timeout=timeout)
        while not r.ok:
            time.sleep(10)
            print(str(datetime.datetime.now())+" : "+str(r.status_code)+" occured. Trying again",file=sys.stderr)
            sys.stderr.flush()
            r = requests.post(request_string,headers=headers,data=data)
        return r.json()
    except Timeout as ex:
        print(str(datetime.datetime.now())+" : Timeout event occured", file=sys.stderr)
        sys.stderr.flush()
        return None


def getChrName(string):
    m=re.search("^NC_0+([^\.0]+)\.*",string)
    if m:
        return m.group(1)
    else:
        return "NA"


def parseSPDI(string):
    m=re.search("^NC_0+([1-9]\d*)\.\d*:(\d+):([ACGTacgt]+):([ACGTacgt]+)",string)
    if m:
        return {"chr":m.group(1),"pos":int(m.group(2))+1,"ref":m.group(3),"alt":m.group(4)}
    else:
        return {"chr":"NA","pos":"NA","ref":"NA","alt":"NA"}


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
server = "http://"+build+".rest.ensembl.org"
if build=="grch38":
    server = "http://rest.ensembl.org"

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
    r=getResponse2(server+ext,headers,list2string(L))
    if r:
        for snprec in r:
            H={}
            id1=snprec["id"][0]
            H[id1]=[]
            spdi=snprec["spdi"]
            for z in spdi:
                m=re.search("^NC_0+",z)
                if m:
                    #print(z)
                    p=parseSPDI(z)
                    H[id1].append(p)
            printHash(H)
    cur_line+=batchsize
