#!/usr/bin/python3

import requests, sys, time
import os
import argparse
import re
from requests.exceptions import Timeout,TooManyRedirects,RequestException
import datetime
import functions

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
    except TooManyRedirects as ex:
        print(str(datetime.datetime.now())+" : TooManyRedirects exception occured", file=sys.stderr)
        sys.stderr.flush()
        return None
    except RequestException as ex:
        print(str(datetime.datetime.now())+" : RequestException occured", file=sys.stderr)
        sys.stderr.flush()
        return None

#----------------------------------------------------------------------------------------------------------------------------------

# default
build="38"
batchsize=200

parser = argparse.ArgumentParser(description="Get chromosome, position, for one rsID or a list of rsIDs",usage="\n\nrs2position.py <INFILE >OUTFILE -b build -s batchsize\nOR\nrs2position.py -r rsID -b build")
parser.add_argument('--build','-b', action="store",help="Genome build: default: 38", default=build)
parser.add_argument('--verbose','-v',default=False,action="store_true",help="verbose output")
parser.add_argument('--size','-s', action="store",help="Batch size: default: 200", default=batchsize)
parser.add_argument('--rs','-r', action="store",help="rsID",required=False)

try:
    args=parser.parse_args()
except:
    sys.exit(0)

if args.build!=None:
    build=args.build

verbose=args.verbose
rsID=args.rs

if args.size!=None:
    batchsize=int(args.size)

ext = "/variant_recoder/homo_sapiens"
server = "https://grch"+build+".rest.ensembl.org"
if build=="38":
    server = "https://rest.ensembl.org"

if sys.stdin.isatty() and not rsID:
    print("No input provided; exit")
    parser.print_help()
    sys.exit(0)

#---------------------------------------------------------------------------------------------------------------------------

IDs=[]

if rsID:
    IDs.append(rsID)
else:
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
            if verbose:
                print("INFO: "+repr(snprec),file=sys.stderr,flush=True)

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
                        H[id0].append(functions.parseSPDI(z))

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
