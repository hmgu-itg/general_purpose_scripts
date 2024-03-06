#!/usr/bin/python3

import argparse
import pandas as pd
import csv
import sys
import os
import re
import logging
import gzip

from itgukbb import utils

def get_log_fname(fname):
    if fname.endswith(".txt"):
        return fname[:-len(".txt")]+".log"
    elif fname.endswith(".txt.gz"):
        return fname[:-len(".txt.gz")]+".log"
    elif fname.endswith(".gz"):
        return fname[:-len(".gz")]+".log"
    else:
        return fname+".log"

# test if second row contains column classes (legacy)
def is_legacy(fname):
    df=pd.read_table(fname,nrows=1,sep="\t",header=0,dtype=str,usecols=["f.eid"],keep_default_na=False)
    return df.iloc[0]["f.eid"]=="NA"

# total number of rows
def get_nrows(fname,legacy=False):
    df=pd.read_table(fname,skiprows=[1] if legacy else None,sep="\t",header=0,dtype=str,usecols=["f.eid"])
    return len(df)

def main():
    verbosity=logging.INFO
    parser=argparse.ArgumentParser(description="This script queries the main phenotype table of a project",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--project','-p',required=True,action='store',help="Project name")
    parser.add_argument('--release','-r',required=True,action='store',help="Project release")
    parser.add_argument('--output','-o',required=False,action='store',help="Output file",default=None)
    parser.add_argument('--chunk','-chunk',required=False,type=int,default=10000,action='store',help="Chunk size; lower it if RAM limit is an issue")
    parser.add_argument('--config','-c',required=False,action='store',help="Config file")
    parser.add_argument('--describe','-d',default=None,metavar="FIELD",required=False,action='store',help="Describe the input data field")
    parser.add_argument('--olink','-olink',required=False,action='store_true',help="Also export OLINK data",default=False)
    parser.add_argument('--list','-l',default=False,required=False,action='store_true',help="Output information about every field in the input project and exit")
    parser.add_argument('--field','-f',metavar="FIELD",required=False,default=[],action='append',help="For a given field, output all its instances; this option can be specified multiple times. Can be a value, or a file containing list of values, or a shell redirection")
    parser.add_argument("--names", "-n",required=False,action='store_true',help="Use field names instead of IDs in output",default=False)
    parser.add_argument("--count-na", "-count-na",required=False,dest="count_na",action='store_true',help="Output NA stats",default=False)
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
    count_na=args.count_na
    to_list=args.list
    d_field=args.describe
    outfname=args.output
    use_names=args.names
    use_olink=args.olink
    chunksize=args.chunk
    output_IDs=False

    logfile=None
    if outfname:
        logfile=get_log_fname(outfname)
    
    if args.verbose is not None:
        if args.verbose=="debug":
            verbosity=logging.DEBUG
        elif args.verbose=="warning":
            verbosity=logging.WARNING
        elif args.verbose=="error":
            verbosity=logging.ERROR

    LOGGER=logging.getLogger("ukbb_select")
    LOGGER.setLevel(verbosity)
    ch1=logging.StreamHandler(sys.stderr)
    ch1.setLevel(verbosity)
    ch2=None
    if logfile:
        ch2=logging.FileHandler(logfile,'w')
    if ch2:
        ch2.setLevel(verbosity)
    formatter=logging.Formatter('%(levelname)s - %(name)s - %(funcName)s -%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    ch1.setFormatter(formatter)
    if ch2:
        ch2.setFormatter(formatter)
    LOGGER.addHandler(ch1)
    if ch2:
        LOGGER.addHandler(ch2)

    logging.getLogger("itgukbb.utils").addHandler(ch1)
    if ch2:
        logging.getLogger("itgukbb.utils").addHandler(ch2)

    if len(args.field)==0 and use_olink==False and d_field is None and to_list==False:
        output_IDs=True
        # LOGGER.error("at least one of the options --olink / --list / --field / --describe must be specified")
        # sys.exit(1)

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
    LOGGER.debug("config file: %s" % config)
    C=utils.readConfig(config)
    if C is None:
        sys.exit(1)

    infile=utils.getProjectFileName(C,project,release,"MAIN")
    if infile is None:
        sys.exit(1)
    
    LOGGER.debug("input file: %s" % infile)
    # key: field, value: dict with keys "Field", "ValueType"
    DICT=utils.readDataDictionary(C["DATA_DICT"])
    olink_fields=list()
    if use_olink:
        olink_fields=[k for k in DICT if DICT[k]["Field"].startswith("OLINK")] # all OLINK fields from data dict

    # analyze the header line
    df=pd.read_table(infile,sep="\t",nrows=1)
    HEADER=dict() # short field name --> list of matching full names: 123 --> [ f.123.1.1, f.123.1.2, ... ]
    for x in df.columns.values.tolist():
        # skip these
        if x=="f.eid" or x=="RELEASE" or x=="CREATED":
            continue
        m=re.match("^f\.(\d+)\.\d+\.\d+$",x)
        if m:
            key=m.group(1)
            HEADER.setdefault(key,[]).append(x)
        else:
            LOGGER.error("column name format error: %s" %x)

    #---------------------------------------------------------------------------------------------------------------------------

    legacy_input=is_legacy(infile)
    LOGGER.info("legacy input: %s" %(str(legacy_input)))

    # just output all sample IDs and exit
    if output_IDs:
        LOGGER.info("output all sample IDs")
        df=pd.read_table(infile,skiprows=[1] if legacy_input else None,sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["f.eid"])
        if outfname:
            df.to_csv(outfname,sep="\t",index=False)
        else:
            print(df.to_csv(sep="\t",index=False),end='')
        sys.exit(0)
        
    #---------------------------------------------------------------------------------------------------------------------------
    
    # describe a field and exit
    if d_field:
        if not d_field in HEADER:
            LOGGER.error("%s is not in input header" % d_field)
            sys.exit(1)
        if not d_field in DICT:
            LOGGER.error("%s is not in data dictionary" % d_field)
            sys.exit(1)
        txt="-" * int(0.85*os.get_terminal_size().columns)
        print(txt.center(os.get_terminal_size().columns))
        print(pd.DataFrame([["Field",d_field],["Field description",DICT[d_field]["Field"]],["Field type",DICT[d_field]["ValueType"]],["Instances",len(HEADER[d_field])]],columns=["A","B"]).to_string(index=False,header=False))
        print(txt.center(os.get_terminal_size().columns))
        if DICT[d_field]["ValueType"]=="Categorical single" or DICT[d_field]["ValueType"]=="Categorical multiple":
            df=pd.read_table(infile,skiprows=[1] if legacy_input else None,sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=HEADER[d_field])
            df=df.apply(pd.Series.value_counts).fillna(0).astype(int)
            df["Count"]=df.agg("sum",axis=1)
            df.index.name="Value"
            df.reset_index(inplace=True)
            print(df.to_string(index=False,columns=["Value","Count"]))
            print(txt.center(os.get_terminal_size().columns))
        sys.exit(0)

    #---------------------------------------------------------------------------------------------------------------------------
    
    # just output info about available fields and exit
    if to_list:
        f=sys.stdout
        if outfname:
            if outfname.endswith(".gz"):
                f=gzip.open(outfname,"wt")
            else:
                f=open(outfname,"w")

        if count_na:
            print("{}\t{}\t{}\t{}\t{}".format("Field","Instances","Description","Type","Missing"),file=f)
            NAcount=dict() # short field name --> # NAs
            for x in HEADER:
                NAcount[x]=0;
            nrows=get_nrows(infile,legacy_input)
            total_chunks=nrows//chunksize
            if nrows%chunksize:
                total_chunks+=1
            LOGGER.info("total rows: %d; total chunks: %d" %(nrows,total_chunks))
            current_chunk=1
            for chunk in pd.read_table(infile,skiprows=[1] if legacy_input else None,sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,chunksize=chunksize):
                LOGGER.info("chunk %d / %d" %(current_chunk,total_chunks))
                for x in HEADER:
                    df=chunk[HEADER[x]]
                    NAcount[x]+=df[df=="NA"].count().sum()
                current_chunk+=1
            LOGGER.info("done")
            for x in HEADER:
                if x in DICT:
                    print("{0}\t{1}\t{2}\t{3}\t{4:.5f}".format(x,len(HEADER[x]),DICT[x]["Field"],DICT[x]["ValueType"],NAcount[x]/(nrows*len(HEADER[x]))),file=f)
                else:
                    print("{}\t{}\t{}\t{}\t{}".format(x,len(HEADER[x]),"NA","NA","NA"),file=f)
                    LOGGER.warning("%s is not in data dictionary" % x)
        else: # not counting NAs
            print("{}\t{}\t{}\t{}".format("Field","Instances","Description","Type"),file=f)
            for x in HEADER:
                if x in DICT:
                    print("{0}\t{1}\t{2}\t{3}".format(x,len(HEADER[x]),DICT[x]["Field"],DICT[x]["ValueType"]),file=f)
                else:
                    print("{}\t{}\t{}\t{}".format(x,len(HEADER[x]),"NA","NA"),file=f)
                    LOGGER.warning("%s is not in data dictionary" % x)            
        sys.exit(0)

    #---------------------------------------------------------------------------------------------------------------------------    
    # general case
    
    LOGGER.info("selecting fields")
    allfields=set()
    for x in args.field:
        allfields.update(utils.readIDList(x))
    allfields.update(olink_fields) # contains short field names
    LOGGER.info("total input fields: %d" %len(allfields))

    # remove fields that are not in header
    tmp=set()
    for f in allfields:
        if not f in HEADER:
            LOGGER.warning("field %s is not in input header" % f)
            tmp.add(f)
    for f in tmp:
        allfields.remove(f)

    to_keep=list() # long field names
    to_keep.append("f.eid")
    for c in allfields:
        for f in HEADER[c]:
            to_keep.append(f)
    to_keep=to_keep[:1]+sorted(to_keep[1:])
    rename_mapper=dict()
    to_keep2=to_keep # long field names with proper name
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
                if s in DICT:
                    s=DICT[s]["Field"]
                else:
                    LOGGER.warn("field %s is not in data dictionary" % s)
                rename_mapper[c]=re.sub(r"^f\.(\d+)(\.\d+\.\d+)$",s+"\\2",c) # f.1234.0.1 --> field_name.0.1
                continue
        to_keep2=list(rename_mapper.values())
        to_keep2=to_keep2[:1]+sorted(to_keep2[1:])
        
    nrows=get_nrows(infile,legacy_input)
    total_chunks=nrows//chunksize
    if nrows%chunksize:
        total_chunks+=1
    LOGGER.info("total rows: %d; total chunks: %d" %(nrows,total_chunks))
        
    current_chunk=1
    with_header=True
    for chunk in pd.read_table(infile,skiprows=[1] if legacy_input else None,sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=to_keep,chunksize=chunksize):
        LOGGER.info("chunk %d / %d" %(current_chunk,total_chunks))
        current_chunk+=1
        chunk.rename(columns=rename_mapper,inplace=True)
        chunk.to_csv(outfname,mode="w" if with_header else "a",sep="\t",columns=to_keep2,index=False,quotechar='"',quoting=csv.QUOTE_NONE,header=with_header)
        with_header=False
    LOGGER.info("done")
    sys.exit(0)

if __name__=="__main__":
    main()
