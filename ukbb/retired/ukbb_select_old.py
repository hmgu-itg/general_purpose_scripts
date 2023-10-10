#!/usr/bin/python3

import argparse
import pandas as pd
import csv
import sys
import os
import re
import logging

# sys.path.append('./python')
from itgukbb import utils

verbosity=logging.INFO

parser=argparse.ArgumentParser(description="This script queries the main phenotype table of a project and outputs information about the provided input fields")
parser.add_argument('--project','-p',required=True,action='store',help="Project name")
parser.add_argument('--release','-r',required=True,action='store',help="Project release")
parser.add_argument('--output','-o',required=False,action='store',help="Output file")
parser.add_argument('--config','-c',required=False,action='store',help="Config file")
parser.add_argument('--describe','-d',metavar="FIELD",required=False,action='store',help="Describe the input data field")
parser.add_argument('--olink','-olink',required=False,action='store_true',help="Also output OLINK data",default=False)
parser.add_argument('--list','-l',required=False,action='store_true',help="Output information about every field in the input project table")
parser.add_argument('--majority','-majority',metavar="FIELD",required=False,action='append',help="For a given field, output the most frequent value, over all instances.  The input field has to have type \"Categorical single\". This option can be specified multiple times")
parser.add_argument('--mean','-mean',metavar="FIELD",required=False,action='append',help="For a given field, output the mean value, over all instances. The input field has to be integer or continuous. This option can be specified multiple times")
parser.add_argument('--min-missing','-min-missing',metavar="FIELD",required=False,action='append',help="For a given field, output the instance with the least number of NAs. This option can be specified multiple times",dest="min_missing")
parser.add_argument('--all','-a',metavar="FIELD",required=False,action='append',help="For a given field, output all its instances. This option can be specified multiple times")
parser.add_argument('--cc','-cc',metavar="FIELD:<COMMA SEPARATED LIST OF CASE VALUES>",required=False,action='append',help="For a given field and a list of \"case\" values, output its case/control (1/0) encoding. The input field has to have type \"Categorical single\" or \"Categorical multiple\". This option can be specified multiple times")
parser.add_argument("--names", "-n",required=False,action='store_true',help="Use field names instead of IDs; default: false")
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
d_field=args.describe
outfname=args.output
use_names=args.names
use_olink=args.olink

if args.verbose is not None:
    if args.verbose=="debug":
        verbosity=logging.DEBUG
    elif args.verbose=="warning":
        verbosity=logging.WARNING
    elif args.verbose=="error":
        verbosity=logging.ERROR

LOGGER=logging.getLogger("ukbb_select")
LOGGER.setLevel(verbosity)
ch=logging.StreamHandler(sys.stderr)
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(funcName)s -%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)

logging.getLogger("itgukbb.utils").addHandler(ch)
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
# key: field, value: dict with keys "Field", "ValueType"
D=utils.readDataDictionary(C["DATA_DICT"])
olink_fields=list()
if use_olink:
    olink_fields=[k for k in D if D[k]["Field"].startswith("OLINK")]

# only the header line
df=pd.read_table(infile,sep="\t",nrows=1)
H=dict() # short field name --> list of matching full names: 123 --> [ f.123.1.1, f.123.1.2, ... ]
for x in df.columns.values.tolist():
    # skip these
    if x=="f.eid" or x=="RELEASE" or x=="CREATED":
        continue
    m=re.match("^f\.(\d+)\.\d+\.\d+$",x)
    if m:
        key=m.group(1)
        H.setdefault(key,[]).append(x)
    else:
        LOGGER.error("column name format error: %s" %x)

#-----------------------------------------------------------------------------------------------------------------------------

# describe a field and exit
if d_field:
    if not d_field in H:
        LOGGER.error("%s is not in input header" % d_field)
        sys.exit(1)
    if not d_field in D:
        LOGGER.error("%s is not in data dictionary" % d_field)
        sys.exit(1)
    txt="-" * int(0.85*os.get_terminal_size().columns)
    print(txt.center(os.get_terminal_size().columns))
    print(pd.DataFrame([["Field",d_field],["Field description",D[d_field]["Field"]],["Field type",D[d_field]["ValueType"]],["Instances",len(H[d_field])]],columns=["A","B"]).to_string(index=False,header=False))
    print(txt.center(os.get_terminal_size().columns))
    if D[d_field]["ValueType"]=="Categorical single" or D[d_field]["ValueType"]=="Categorical multiple":
        df=pd.read_table(infile,skiprows=[1],sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=H[d_field])
        df=df.apply(pd.Series.value_counts).fillna(0).astype(int)
        df["Count"]=df.agg("sum",axis=1)
        df.index.name="Value"
        df.reset_index(inplace=True)
        print(df.to_string(index=False,columns=["Value","Count"]))
        print(txt.center(os.get_terminal_size().columns))
    sys.exit(0)

# just output info about available fields and exit
if to_list:
    f=sys.stdout
    if outfname:
        f=open(outfname,"w")
    print("{}\t{}\t{}\t{}".format("Field","Instances","Description","Type"),file=f)
    for x in H:
        if x in D:
            print("{}\t{}\t{}\t{}".format(x,len(H[x]),D[x]["Field"],D[x]["ValueType"]),file=f)
        else:
            print("{}\t{}\t{}\t{}".format(x,len(H[x]),"NA","NA"),file=f)
            LOGGER.warning("%s is not in data dictionary" % x)
    sys.exit(0)
else:
    LOGGER.info("selecting fields")
    ccfields=set() # set of f:v1,v2,v3 ... strings
    if not args.cc is None:
        for x in args.cc:
            ccfields.add(x)
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
    allfields.update(olink_fields)

    # check field availability
    tmp=set()
    for f in ccfields:
        x=f.split(":")[0]
        if x not in H:
            LOGGER.warning("field %s is not in input header" % x)
            tmp.add(f)
    for f in tmp:
        ccfields.remove(f)
        
    tmp.clear()
    for s in [majfields,meanfields,minmissfields,allfields]:
        for f in s:
            if not f in H:
                LOGGER.warning("field %s is not in input header" % f)
                tmp.add(f)
        for f in tmp:
            s.remove(f)
        tmp.clear()

    # check field types
    tmp.clear()
    for f in ccfields:
        x=f.split(":")[0]
        t=D[x]["ValueType"]
        if t!="Categorical single" and t!="Categorical multiple":
            LOGGER.warning("case/control field %s has type %s; only \"Categorical single\" or \"Categorical multiple\" type is allowed for case/control" % (x,t))
            tmp.add(f)
    for x in tmp:
        ccfields.remove(x)
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

#-----------------------------------------------------------------------------------------------------------------------------

    to_keep=set()
    to_keep.add("f.eid")
    for f in ccfields:
        z=f.split(":")[0]
        for x in H[z]:
            to_keep.add(x)
    for s in [majfields,meanfields,minmissfields,allfields]:
        for f in s:
            for x in H[f]:
                to_keep.add(x)
    # LOGGER.info("reading columns %s" % ", ".join(str(e) for e in to_keep))
    # skip second row
    df=pd.read_table(infile,skiprows=[1],sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=to_keep)
    
    to_keep=list()
    to_keep.append("f.eid")
    for c in allfields:
        for f in H[c]:
            to_keep.append(f)
    for c in ccfields:
        z=c.split(":")
        LOGGER.info("creating CC column for %s and value(s) %s" %(z[0],z[1]))
        new_colname="cc-"+z[0]+"-"+"_".join(str(e) for e in z[1].split(","))
        df=utils.addSummaryColumn(df,H[z[0]],new_colname,z[1].split(","),"cc")
        if not new_colname in to_keep:
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
    #LOGGER.info("writing columns %s" % ", ".join(str(e) for e in to_keep))
    rename_mapper=dict()
    if use_names:
        for c in to_keep:
            if c=="f.eid":
                rename_mapper[c]="eid"
                continue
            m=re.match("^f\.(\d+)\.\d+\.\d+$",c)
            if m:
                s=m.group(1)
                if s in D:
                    s=D[s]["Field"]
                else:
                    LOGGER.warn("Field %s is not in data dictionary" % s)
                rename_mapper[c]=re.sub(r"^f\.(\d+)(\.\d+\.\d+)$",s+"\\2",c)
                continue
            m=re.match("^cc-(\d+)-.*$",c)
            if m:
                s=m.group(1)
                if s in D:
                    s=D[s]["Field"]
                else:
                    LOGGER.warn("Field %s is not in data dictionary" % s)
                rename_mapper[c]=re.sub(r"^cc-(\d+)-(.*)$","cc-"+s+"-\\2",c)
                continue
            m=re.match("^majority-(\d+)$",c)
            if m:
                s=m.group(1)
                if s in D:
                    s=D[s]["Field"]
                else:
                    LOGGER.warn("Field %s is not in data dictionary" % s)
                rename_mapper[c]="majority-"+s
                continue
            m=re.match("^mean-(\d+)$",c)
            if m:
                s=m.group(1)
                if s in D:
                    s=D[s]["Field"]
                else:
                    LOGGER.warn("Field %s is not in data dictionary" % s)
                rename_mapper[c]="mean-"+s
                continue
            m=re.match("^min_missing-(\d+)$",c)
            if m:
                s=m.group(1)
                if s in D:
                    s=D[s]["Field"]
                else:
                    LOGGER.warn("Field %s is not in data dictionary" % s)
                rename_mapper[c]="min_missing-"+s
                continue
        to_keep=list(rename_mapper.values())
    else:
        # just rename f.eid --> eid
        rename_mapper["f.eid"]="eid"
        to_keep=["eid" if x=="f.eid" else x for x in to_keep]
    LOGGER.debug(rename_mapper)
    df.rename(columns=rename_mapper,inplace=True)
    if outfname:
        df.to_csv(outfname,sep="\t",columns=to_keep,index=False,quotechar='"',quoting=csv.QUOTE_NONE)
    else:
        print(df.to_string(index=False,columns=to_keep))        
