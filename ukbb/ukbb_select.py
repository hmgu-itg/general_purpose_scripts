#!/usr/bin/python3

import argparse
import pandas as pd
import csv
import sys
import os
from itgukbb import utils
import re
import logging

verbosity=logging.INFO

parser=argparse.ArgumentParser()
parser.add_argument('--project','-p',required=True,action='store',help="Project name")
parser.add_argument('--release','-r',required=True,action='store',help="Release")
parser.add_argument('--output','-o',required=True,action='store',help="Output prefix")
parser.add_argument('--config','-c',required=False,action='store',help="Config file")
parser.add_argument('--list','-l',required=False,action='store_true',help="Output column information")
parser.add_argument('--majority','-majority',required=False,action='append',help="Output most frequent value, over all instances")
parser.add_argument('--mean','-mean',required=False,action='append',help="Output mean, over all instances")
parser.add_argument('--min-missing','-min-missing',required=False,action='append',help="Output instance with least NAs",dest="min_missing")
parser.add_argument('--all','-a',required=False,action='append',help="Output all instances")
parser.add_argument('--cc','-cc',required=False,action='append',help="Recode as case/control (1/0)")
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
to_list=args.list
out_prefix=args.output
logF=out_prefix+".log"

if args.verbose is not None:
    if args.verbose=="debug":
        verbosity=logging.DEBUG
    elif args.verbose=="warning":
        verbosity=logging.WARNING
    elif args.verbose=="error":
        verbosity=logging.ERROR

LOGGER=logging.getLogger("ukbb_select")
LOGGER.setLevel(verbosity)
ch=logging.FileHandler(logF,'w')
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(funcName)s -%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)
LOGGER.addHandler(logging.StreamHandler(sys.stdout))

logging.getLogger("itgukbb.utils").addHandler(ch)
logging.getLogger("itgukbb.utils").addHandler(logging.StreamHandler(sys.stdout))
logging.getLogger("itgukbb.utils").setLevel(verbosity)

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

infile=utils.getProjectFileName(C,project,release,"MAIN")
if infile is None:
    sys.exit(1)
    
LOGGER.info("input file: %s" % infile)
# key: field
D=utils.readDataDictionary(C["DATA_DICT"])

df=pd.read_table(infile,sep="\t",nrows=1)
H=dict() # short field name --> list of matching full names
for x in df.columns.values.tolist():
    if x=="f.eid" or x=="RELEASE" or x=="CREATED":
        continue
    m=re.match("^f\.(\d+)\.\d+\.\d+$",x)
    if m:
        key=m.group(1)
        H.setdefault(key,[]).append(x)
    else:
        LOGGER.error("column name format error: %s" %x)

#-----------------------------------------------------------------------------------------------------------------------------

if to_list:
    LOGGER.info("output field info")
    with open(out_prefix+".txt","w") as f:
        print("{}\t{}\t{}\t{}".format("Field","Instances","Description","Type"),file=f)
        for x in H:
            if x in D:
                print("{}\t{}\t{}\t{}".format(x,len(H[x]),D[x]["Field"],D[x]["ValueType"]),file=f)
            else:
                print("{}\t{}\t{}\t{}".format(x,len(H[x]),"NA","NA"),file=f)
                LOGGER.warning("%s is not in data dictionary" % x)
else:
    LOGGER.info("selecting fields")
    ccfields=dict() # field --> set of "case" values
    if not args.cc is None:
        for x in args.cc:
            a=x.split(":",1)
            if len(a)!=2:
                LOGGER.error("wrong --cc option value: %s" % x)
                continue
            b=a[1].split(",")
            for z in b:
                ccfields.setdefault(a[0],set()).add(z)
    majfields=set()
    if not args.majority is None:
        for x in args.majority:
            majfields.add(x)
    meanfields=set()
    if not args.mean is None:
        for x in args.mean:
            meanfields.add(x)
    minmissfields=set()
    if not args.min_missing is None:
        for x in args.min_missing:
            minmissfields.add(x)
    allfields=set()
    if not args.all is None:
        for x in args.all:
            allfields.add(x)

    # check field availability
    tmp=set()
    for s in [ccfields,majfields,meanfields,minmissfields,allfields]:
        for f in s:
            if not f in H:
                LOGGER.warning("field %s is not in input header" % f)
                tmp.add(f)
        for f in tmp:
            if s==ccfields:
                del s[f]
            else:
                s.remove(f)
        tmp.clear()

    # check field types
    tmp.clear()
    for f in ccfields:
        t=D[f]["ValueType"]
        if t!="Categorical single" and t!="Categorical multiple":
            LOGGER.warning("case/control field %s has type %s; only \"Categorical single\" or \"Categorical multiple\" type is allowed for case/control" % (f,t))
            tmp.add(f)
    for x in tmp:
        del ccfields[x]
    tmp.clear()

    for f in majfields:
        t=D[f]["ValueType"]
        if t!="Categorical single":
            LOGGER.warning("majority field %s has type %s; only \"Categorical single\" type is allowed for majority" % (f,t))
            tmp.add(f)
    for x in tmp:
        del majfields[x]
    tmp.clear()

    for f in meanfields:
        t=D[f]["ValueType"]
        if t!="Continuous" and t!="Integer":
            LOGGER.warning("mean field %s has type %s; only \"Continuous\" and \"Integer\" types are allowed for mean" % (f,t))
            tmp.add(f)
    for x in tmp:
        del meanfields[x]
    tmp.clear()

    to_keep=set()
    to_keep.add("f.eid")
    for s in [ccfields,majfields,meanfields,minmissfields,allfields]:
        for f in s:
            for x in H[f]:
                to_keep.add(x)

#-----------------------------------------------------------------------------------------------------------------------------

    LOGGER.info("reading columns %s" % ", ".join(str(e) for e in to_keep))
    df=pd.read_table(infile,skiprows=[1],sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=to_keep)
    to_keep=list()
    to_keep.append("f.eid")
    for c in allfields:
        for f in H[c]:
            to_keep.append(f)
    for c in ccfields:
        LOGGER.info("creating CC column for %s" %c)
        s="_".join(str(e) for e in ccfields[c])
        new_colname="cc-"+c+"-"+s
        df=utils.addSummaryColumn(df,H[c],new_colname,list(ccfields[c]),"cc")
        to_keep.append(new_colname)
    for c in majfields:
        LOGGER.info("creating majority column for %s" %c)
        new_colname="majority-"+c
        df=utils.addSummaryColumn(df,H[c],new_colname,None,"majority")
        to_keep.append(new_colname)
    for c in meanfields:
        LOGGER.info("creating mean column for %s" %c)
        new_colname="mean-"+c
        df=utils.addSummaryColumn(df,H[c],new_colname,None,"mean")
        to_keep.append(new_colname)
    for c in minmissfields:
        LOGGER.info("creating min missing column for %s" %c)
        new_colname="min_missing-"+c
        df=utils.addSummaryColumn(df,H[c],new_colname,None,"minmissing")
        to_keep.append(new_colname)
    LOGGER.info("writing columns %s" % ", ".join(str(e) for e in to_keep))
    df.to_csv(out_prefix+".txt",sep="\t",columns=to_keep,index=False,quotechar='"',quoting=csv.QUOTE_NONE)
