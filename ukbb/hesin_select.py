#!/usr/bin/python3

import argparse
import pandas as pd
import csv
import sys
import tarfile
import logging
import os
from itgukbb import utils

verbosity=logging.INFO

parser=argparse.ArgumentParser()
parser.add_argument('--project','-p',required=True,action='store',help="Project name")
parser.add_argument('--release','-r',required=True,action='store',help="Release")
parser.add_argument('-icd9','--icd9',required=False,action='store',help="List of ICD9 codes")
parser.add_argument('-icd10','--icd10',required=False,action='store',help="List of ICD10 codes")
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

icd9codes=list()
icd10codes=list()

if not args.icd9 is None:
    with open(args.icd9,"r") as f:
        icd9codes=f.read().splitlines()

if not args.icd10 is None:
    with open(args.icd10,"r") as f:
        icd10codes=f.read().splitlines()

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

L=list()
with tarfile.open(infile,"r:*") as tar:
    df=pd.read_table(tar.extractfile("hesin_diag.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","diag_icd9","diag_icd10"])
    if icd10codes:
        L=df.loc[df["diag_icd10"].isin(icd10codes)]["eid"].unique().tolist()
    if icd9codes:
        L.extend(x for x in df.loc[df["diag_icd9"].isin(icd9codes)]["eid"].unique().tolist() if not x in L)
    LOGGER.info("output %d ID(s)" %len(L))
    with open(outfname,"w") as f:
        if len(L):
            print("\n".join(L),file=f)
