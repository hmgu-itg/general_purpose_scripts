#!/usr/bin/python3

import argparse
import pandas as pd
import csv
import sys
import tarfile
import logging
import os
import re
from itgukbb import utils
import functools as ft

def tramsformExpr(string):
    return re.sub(r"(\w+)",lambda x:"(\""+x.group(1)+"\" in L)" if (x.group(1).lower()!="and" and x.group(1).lower()!="or") else x.group(1),string)

def testF(L,L2):
    return ft.reduce(lambda a,b:a or b,list(map(lambda x:eval(x),L2)))

verbosity=logging.INFO

parser=argparse.ArgumentParser()
parser.add_argument('--project','-p',required=True,action='store',help="Project name")
parser.add_argument('--release','-r',required=True,action='store',help="Release")
parser.add_argument('-icd9','--icd9',required=False,action='store',help="List of ICD9 codes")
parser.add_argument('-icd10','--icd10',required=False,action='store',help="List of ICD10 codes")
parser.add_argument('-oper3','--oper3',required=False,action='store',help="List of OPCS3 codes")
parser.add_argument('-oper4','--oper4',required=False,action='store',help="List of OPCS4 codes")
parser.add_argument('-o','--output',required=True,action='store',help="Output file")
parser.add_argument('--config','-c',required=False,action='store',help="Config file")
parser.add_argument("--verbose", "-v", help="Verbosity level; default: info",required=False,choices=("debug","info","warning","error"),default="info")

if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

try:
    args=parser.parse_args()
except:
    sys.exit(1)

if args.icd9 is None and args.icd10 is None:
    print("ERROR: no ICD codes given",file=sys.stderr)
    sys.exit(1)
    
project=args.project
release=args.release
outfname=args.output

if args.verbose is not None:
    if args.verbose=="debug":
        verbosity=logging.DEBUG
    elif args.verbose=="warning":
        verbosity=logging.WARNING
    elif args.verbose=="error":
        verbosity=logging.ERROR

LOGGER=logging.getLogger("hesin_select")
LOGGER.setLevel(verbosity)
ch=logging.StreamHandler(sys.stderr)
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(funcName)s -%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("itgukbb.utils").addHandler(ch)
logging.getLogger("itgukbb.utils").setLevel(verbosity)

icd9codes=set()
icd10codes=set()
opcs3codes=set()
opcs4codes=set()

if not args.icd9 is None:
    with open(args.icd9,"r") as f:
        for l in f.read().splitlines():
            if not re.match("^\s*$",l) and not re.match("^#.*",l):
                icd9codes.add(l)

if not args.icd10 is None:
    with open(args.icd10,"r") as f:
        for l in f.read().splitlines():
            if not re.match("^\s*$",l) and not re.match("^#.*",l):
                icd10codes.add(l)

if not args.opcs3 is None:
    with open(args.opcs3,"r") as f:
        for l in f.read().splitlines():
            if not re.match("^\s*$",l) and not re.match("^#.*",l):
                opcs3codes.add(l)

if not args.opcs4 is None:
    with open(args.opcs4,"r") as f:
        for l in f.read().splitlines():
            if not re.match("^\s*$",l) and not re.match("^#.*",l):
                opcs4codes.add(l)

# TODO: at least one list should be provided
        
#-----------------------------------------------------------------------------------------------------------------------------

if args.config is None:
    config=os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),"config.txt")
else:
    config=args.config

for arg in vars(args):
    LOGGER.info("INPUT OPTIONS: %s : %s" % (arg, getattr(args, arg)))
LOGGER.info("")
LOGGER.info("config file: %s" % config)
C=utils.readConfig(config)
if C is None:
    sys.exit(1)

infile=utils.getProjectFileName(C,project,release,"HESIN")
if infile is None:
    sys.exit(1)

LOGGER.info("input file: %s" % infile)
LOGGER.info("output file: %s" % outfname)

#-----------------------------------------------------------------------------------------------------------------------------

Licd9=set()
Licd10=set()
Loper3=set()
Loper4=set()
L=list()
with tarfile.open(infile,"r:*") as tar:
    df_diag=pd.read_table(tar.extractfile("hesin_diag.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","ins_index","diag_icd9","diag_icd10"])
    df_oper=pd.read_table(tar.extractfile("hesin_oper.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","ins_index","oper3","oper4"])
    df_diag2=df_diag.groupby(["eid","ins_index"],as_index=False).agg({"diag_icd9":lambda x:list(x),"diag_icd10":lambda x:list(x)})
    df_oper2=df_oper.groupby(["eid","ins_index"],as_index=False).agg({"oper3":lambda x:list(x),"oper4":lambda x:list(x)})
    if icd10codes:
        Licd10=set(df_diag2[df_diag2.apply(lambda x:testF(x["diag_icd10"],[transformExpr(z) for z in icd10codes]))==True,axis=1]["eid"].tolist())
    if icd9codes:
        Licd9=set(df_diag2[df_diag2.apply(lambda x:testF(x["diag_icd9"],[transformExpr(z) for z in icd9codes]))==True,axis=1]["eid"].tolist())
    if opcs3codes:
        Loper3=set(df_oper2[df_oper2.apply(lambda x:testF(x["oper3"],[transformExpr(z) for z in opcs3codes]))==True,axis=1]["eid"].tolist())
    if opcs4codes:
        Loper4=set(df_oper2[df_oper2.apply(lambda x:testF(x["oper4"],[transformExpr(z) for z in opcs4codes]))==True,axis=1]["eid"].tolist())
    L=list(Licd10.union(Licd9,Loper3,Loper4))
    LOGGER.info("output %d ID(s)" %len(L))
    with open(outfname,"w") as f:
        if len(L):
            print("\n".join(L),file=f)
