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

def transformExpr(string):
    return re.sub(r"(\w+)",lambda x:"(\""+x.group(1)+"\" in L)" if (x.group(1).lower()!="and" and x.group(1).lower()!="or") else x.group(1).lower(),string)

def testF(L1,L2):
    return ft.reduce(lambda a,b:a or b,list(map(lambda x:eval(x,{},{"L":L1}),L2)))

def testF2(L1,L2):
    for x in L1:
        if ft.reduce(lambda a,b:a or b,list(map(lambda c:eval(c,{},{"L":x}),L2))):
            return True
    return False

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

#----------------------------------------------------------------------------------------------------------------------------------

if args.icd9 is None and args.icd10 is None and args.oper3 is None and args.oper4 is None:
    LOGGER.error("no ICD/OPCS codes given")
    sys.exit(1)

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

if not args.oper3 is None:
    with open(args.oper3,"r") as f:
        for l in f.read().splitlines():
            if not re.match("^\s*$",l) and not re.match("^#.*",l):
                opcs3codes.add(l)

if not args.oper4 is None:
    with open(args.oper4,"r") as f:
        for l in f.read().splitlines():
            if not re.match("^\s*$",l) and not re.match("^#.*",l):
                opcs4codes.add(l)

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

# assuming ICD lists don't contain logical expressions

Licd9=set()
Licd10=set()
Loper3=set()
Loper4=set()
L=list()
with tarfile.open(infile,"r:*") as tar:
    if icd9codes or icd10codes:
        df_diag=pd.read_table(tar.extractfile("hesin_diag.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","ins_index","diag_icd9","diag_icd10"])
        if icd10codes:
            Licd10=set(df_diag.loc[df_diag["diag_icd10"].isin(icd10codes)]["eid"].tolist())
            LOGGER.info("%d IDs match ICD10 codes" % len(Licd10))
        if icd9codes:
            Licd9=set(df_diag.loc[df_diag["diag_icd9"].isin(icd9codes)]["eid"].tolist())
            LOGGER.info("%d IDs match ICD9 codes" % len(Licd9))

    if opcs3codes or opcs4codes:
        df_oper=pd.read_table(tar.extractfile("hesin_oper.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","ins_index","oper3","oper4"])
        df_oper2=df_oper.groupby(["eid","ins_index"],as_index=False).agg({"oper3":lambda x:list(x),"oper4":lambda x:list(x)})
        df_oper2=df_oper2[["eid","oper3","oper4"]]
        df_oper2=df_oper2.groupby(["eid"],as_index=False).agg({"oper3":lambda x:list(x),"oper4":lambda x:list(x)})
        L3=[transformExpr(z) for z in opcs3codes]
        L4=[transformExpr(z) for z in opcs4codes]
        if opcs3codes:
            # Loper3=set(df_oper2[df_oper2.apply(lambda x:testF(x["oper3"],L3)==True,axis=1)]["eid"].tolist())
            Loper3=set(df_oper2[df_oper2.apply(lambda x:testF2(x["oper3"],L3)==True,axis=1)]["eid"].tolist())
            LOGGER.info("%d IDs match OPCS3 codes" % len(Loper3))
        if opcs4codes:
            # Loper4=set(df_oper2[df_oper2.apply(lambda x:testF(x["oper4"],L4)==True,axis=1)]["eid"].tolist())
            Loper4=set(df_oper2[df_oper2.apply(lambda x:testF2(x["oper4"],L4)==True,axis=1)]["eid"].tolist())
            LOGGER.info("%d IDs match OPCS4 codes" % len(Loper4))
        
    L=list(Licd10.union(Licd9,Loper3,Loper4))
    LOGGER.info("output %d ID(s)" %len(L))
    with open(outfname,"w") as f:
        if len(L):
            print("\n".join(L),file=f)
