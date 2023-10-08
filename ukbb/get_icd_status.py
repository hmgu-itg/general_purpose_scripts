#!/usr/bin/python3

import argparse
import pandas as pd
import csv
import sys
import tarfile
import logging
import os
import re
import stat
from itgukbb import utils

# given patient ID(s) and ICD code, get Y/N if a patient was ever deagnosed with this ICD

# string can be either ID, or filename, or name of a pipe
def readIDList(string):
    if string is None:
        return None
    if os.path.isfile(string) or stat.S_ISFIFO(os.stat(string).st_mode):
        df=pd.read_csv(string,header=None,usecols=[0])
        return list(set(df[0].astype(str).tolist())) # keep unique IDs
    else:
        return [string]

def main():
    verbosity=logging.INFO

    parser=argparse.ArgumentParser(description="Given patient IDs and ICD code, get patient's 0/1 status",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    rgroup=parser.add_argument_group("required arguments")
    rgroup.add_argument('--project','-p',required=True,action='store',help="Project name")
    rgroup.add_argument('--release','-r',required=True,action='store',help="Project release")
    group=rgroup.add_mutually_exclusive_group(required=True)
    group.add_argument('--icd10','-icd10',required=False,action='store',help="ICD10 code")
    group.add_argument('--icd9','-icd9',required=False,action='store',help="ICD9 code")
    parser.add_argument('--id','-id',required=False,action='store',help="Patient ID or a file with patient IDs")
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
    id_list=readIDList(args.id)
    icd10=args.icd10
    icd9=args.icd9

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
        df=df.groupby(["eid"],as_index=False).agg({"icd":lambda x:list(x)})
        df[icd]=df.apply(lambda x: 1 if icd in x["icd"] else 0,axis=1)
        print(df[["eid",icd]].rename(columns={"eid":"ID"}).to_csv(sep="\t",index=False),end='') 
        
if __name__=="__main__":
    main()
