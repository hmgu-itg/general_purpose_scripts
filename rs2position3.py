#!/usr/bin/python3

import requests, sys, time
import os
import argparse
import re
import time
from requests.exceptions import Timeout,TooManyRedirects,RequestException
import datetime

headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

def list2string(snps):
    return "{\"ids\":["+",".join(list(map(lambda x:"\""+x+"\"",snps)))+"]}"

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
    except requests.exceptions.TooManyRedirects:
        print(str(datetime.datetime.now())+" : TooManyRedirects exception occured", file=sys.stderr)
        sys.stderr.flush()
        return None
    except requests.exceptions.RequestException as e:
        print(str(datetime.datetime.now())+" : RequestException occured", file=sys.stderr)
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

#----------------------------------------------------------------------------------------------------------------------------------

build="grch38"
batchsize=200

parser = argparse.ArgumentParser(description="Get chromosome, position, REF and ALT alleles for a list of rsIDs\nINPUT: STDIN\nOUTPUT: STDOUT\n./rs2position3.py <INFILE >OUTFILE")
parser.add_argument('--build','-b', action="store",help="Genome build: default: grch38", default="grch38")
parser.add_argument('--size','-s', action="store",help="Batch size: default: 200", default=200)

try:
    args=parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

if args.build!=None:
    build=args.build

if args.size!=None:
    batchsize=int(args.size)

ext = "/variant_recoder/homo_sapiens"
server = "https://"+build+".rest.ensembl.org"
if build=="grch38":
    server = "https://rest.ensembl.org"

#---------------------------------------------------------------------------------------------------------------------------

IDs=sys.stdin.readlines()
total_lines=len(IDs)

cur_line=0
while cur_line<total_lines:
    L=[]
    for i in range(cur_line,min(cur_line+batchsize,total_lines)):
        x=IDs[i].rstrip(os.linesep)
        if x not in L:
            L.append(x)    

    r=getResponse2(server+ext,headers,list2string(L),max_attempts=5)
    if r:
        for snprec in r:
            H={}
            id0=snprec["input"] # original ID
            L.remove(id0)
            id1=snprec["id"][0]
            if id1!=id0:
                print("WARNING: INPUT ID="+id0,"RETRIEVED ID="+id1,file=sys.stderr,flush=True)

            H[id0]=[]
            if "spdi" in snprec:
                spdi=snprec["spdi"]
                for z in spdi:
                    m=re.search("^NC_0+",z)
                    if m:
                        H[id0].append(parseSPDI(z))

            s=H[id0]
            positions=set(x["chr"]+":"+str(x["pos"]) for x in s)
            if len(positions)>1:
                print("ERROR: more than one position for "+id0,file=sys.stderr,flush=True)
            elif len(positions)==0:
                print("ERROR: no position for "+id0,file=sys.stderr,flush=True)
                print(id0,"NA","NA",sep='\t',file=sys.stdout,flush=True)    
            else:
                L1=positions.pop().rsplit(":")
                print(id0,L1[0],L1[1],sep='\t',file=sys.stdout,flush=True)

    for x in L:
        print(x,"NA","NA",sep='\t',file=sys.stdout,flush=True)    
                    
    cur_line+=batchsize
