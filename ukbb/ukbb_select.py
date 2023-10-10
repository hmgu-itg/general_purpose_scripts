#!/usr/bin/python3

import argparse
import pandas as pd
import csv
import sys
import os
import re
import logging

from itgukbb import utils

# test if second row contains column classes (legacy)
def is_legacy(fname):
    df=pd.read_table(fname,nrows=1,sep="\t",header=0,dtype=str,usecols=["f.eid"],keep_default_na=False)
    return df.iloc[0]["f.eid"]=="NA"

def get_nrows(fname,legacy=False):
    df=pd.read_table(fname,skiprows=[1] if legacy else None,sep="\t",header=0,dtype=str,usecols=["f.eid"])
    return len(df)

verbosity=logging.INFO

parser=argparse.ArgumentParser(description="This script queries the main phenotype table of a project",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--project','-p',required=True,action='store',help="Project name")
parser.add_argument('--release','-r',required=True,action='store',help="Project release")
parser.add_argument('--output','-o',required=False,action='store',help="Output file",default=None)
parser.add_argument('--chunk','-chunk',required=False,type=int,default=10000,action='store',help="Chunk size; lower it if RAM limit is an issue")
parser.add_argument('--config','-c',required=False,action='store',help="Config file")
parser.add_argument('--describe','-d',default=None,metavar="FIELD",required=False,action='store',help="Describe the input data field")
parser.add_argument('--olink','-olink',required=False,action='store_true',help="Also export OLINK data",default=False)
parser.add_argument('--list','-l',default=False,required=False,action='store_true',help="Output information about every field in the input project")
parser.add_argument('--field','-f',metavar="FIELD",required=False,default=[],action='append',help="For a given field, output all its instances. This option can be specified multiple times")
parser.add_argument("--names", "-n",required=False,action='store_true',help="Use field names instead of IDs in output",default=False)
parser.add_argument("--verbose", "-v", help="Verbosity level",required=False,choices=("debug","info","warning","error"),default="info")

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
chunksize=args.chunk

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

if len(args.field)==0 and use_olink==False and d_field is None and to_list==False:
    LOGGER.error("at least one of the options --olink / --list / --field / --describe must be specified")
    sys.exit(1)

if (len(args.field)!=0 or use_olink) and outfname is None:
    LOGGER.error(" --output is not specified")
    sys.exit(1)
    
#-------------------------------------------------------------------------------------------------------------------------

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

# analyze the header line
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

#---------------------------------------------------------------------------------------------------------------------------

legacy_input=is_legacy(infile)
LOGGER.info("legacy input: %s" %(str(legacy_input)))

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
        df=pd.read_table(infile,skiprows=[1] if legacy_input else None,sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=H[d_field])
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
    allfields=set()
    for x in args.field:
        allfields.add(x)
    allfields.update(olink_fields)

    # check field availability
    tmp=set()
    for f in allfields:
        if not f in H:
            LOGGER.warning("field %s is not in input header" % f)
            tmp.add(f)
    for f in tmp:
        allfields.remove(f)

#---------------------------------------------------------------------------------------------------------------------------
    
    to_keep=list()
    to_keep.append("f.eid")
    for c in allfields:
        for f in H[c]:
            to_keep.append(f)
    rename_mapper=dict()
    to_keep2=to_keep
    if use_names:
        to_keep2=list()
        to_keep2.append("eid")
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
        to_keep2=list(rename_mapper.values())
    else:
        to_keep2=to_keep
    c=1
    with_header=True
    nrows=get_nrows(infile,legacy_input)
    LOGGER.debug("Rows: %d" %(nrows))
    tc=nrows//chunksize
    if nrows%chunksize:
        tc+=1
    for chunk in pd.read_table(infile,skiprows=[1] if legacy_input else None,sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=to_keep,chunksize=chunksize):
        LOGGER.info("chunk %d / %d" %(c,tc))
        c+=1
        chunk.rename(columns=rename_mapper,inplace=True)
        if with_header:
            out_mode="w"
        else:
            out_mode="a"
        chunk.to_csv(outfname,mode=out_mode,sep="\t",columns=to_keep2,index=False,quotechar='"',quoting=csv.QUOTE_NONE,header=with_header)
        with_header=False
    LOGGER.info("Done")
