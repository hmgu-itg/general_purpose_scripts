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

# given patient ID(s) and ICD code, get Y/N if a patient was ever deagnosed with this ICD

def main():
    verbosity=logging.INFO

    parser=argparse.ArgumentParser(description="Given patient IDs and ICD code(s), get patient's 0/1 status. Status 1 means a patient was at some time point diagnosed with at least one of provided ICD codes, otherwise status is 0",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    rgroup=parser.add_argument_group("required arguments")
    rgroup.add_argument('--project','-p',required=True,action='store',help="Project name")
    rgroup.add_argument('--release','-r',required=True,action='store',help="Project release")
    group=rgroup.add_mutually_exclusive_group(required=True)
    group.add_argument('--icd10','-icd10',required=False,action='store',help="ICD10 code or a file with ICD10 code(s)")
    group.add_argument('--icd9','-icd9',required=False,action='store',help="ICD9 code or a file with ICD9 code(s)")
    parser.add_argument('--id','-id',required=False,action='store',help="Patient ID or a file with patient IDs")
    parser.add_argument('--config','-c',required=False,action='store',help="Config file")
    parser.add_argument('--output','-o',required=False,action='store',help="Output file")
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
    id_list=utils.readIDList(args.id)
    icd10=utils.readIDList(args.icd10)
    icd9=utils.readIDList(args.icd9)

    if args.verbose is not None:
        if args.verbose=="debug":
            verbosity=logging.DEBUG
        elif args.verbose=="warning":
            verbosity=logging.WARNING
        elif args.verbose=="error":
            verbosity=logging.ERROR

    LOGGER=logging.getLogger("get_icd_status")
    LOGGER.setLevel(verbosity)
    ch=logging.StreamHandler(sys.stderr)
    ch.setLevel(verbosity)
    formatter=logging.Formatter('%(levelname)s - %(name)s - %(funcName)s -%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    ch.setFormatter(formatter)
    LOGGER.addHandler(ch)

    logging.getLogger("itgukbb.utils").addHandler(ch)
    logging.getLogger("itgukbb.utils").setLevel(verbosity)

#----------------------------------------------------------------------------------------------------------------------------------

    if args.config is None:
        config=os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),"config.txt")
    else:
        config=args.config

    for arg in vars(args):
        LOGGER.info("INPUT OPTIONS: %s : %s" % (arg, getattr(args, arg)))

    LOGGER.info("")
    LOGGER.info("config file: %s" % config)
    CONFIG=utils.readConfig(config)
    if CONFIG is None:
        sys.exit(1)
    infile=utils.getProjectFileName(CONFIG,project,release,"HESIN")
    if infile is None:
        sys.exit(1)

    LOGGER.info("input file: %s" % infile)
    if id_list is None or len(id_list)==0:
        LOGGER.info("using all samples in input file")
    LOGGER.info("")

    icd=icd9
    icd_col="diag_icd9"
    if icd9 is None:
        icd=icd10
        icd_col="diag_icd10"
        
#-----------------------------------------------------------------------------------------------------------------------------

    with tarfile.open(infile,"r:*") as tar:
        df=pd.read_table(tar.extractfile("hesin_diag.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid",icd_col])
    if not(id_list is None) and len(id_list)!=0:
        df=df[df["eid"].isin(id_list)]
    if len(df)==0:
        LOGGER.info("no records left after selecting sample IDs from the input list")
        sys.exit(0)
    df=df.groupby(["eid"],as_index=False).agg({icd_col:lambda x:list(x)})
    df["status"]=df.apply(lambda x: 0 if set(icd).isdisjoint(set(x[icd_col])) else 1,axis=1)
    if args.output:
        df[["eid","status"]].rename(columns={"eid":"ID"}).to_csv(args.output,sep="\t",index=False)
    else:
        print(df[["eid","status"]].rename(columns={"eid":"ID"}).to_csv(sep="\t",index=False),end='')        
    if not(id_list is None) and len(id_list)!=0:
        # IDs in input list but not in the previous output
        not_found_ids=set(id_list)-set(df["eid"])
        if len(not_found_ids)!=0:
            if args.output:
                pd.DataFrame({"A":list(not_found_ids),"B":"NA"}).to_csv(args.output,header=False,index=False,sep="\t")
            else:
                print(pd.DataFrame({"A":list(not_found_ids),"B":"NA"}).to_csv(header=False,index=False,sep="\t"),end='')
                
if __name__=="__main__":
    main()
