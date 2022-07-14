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

# given patient ID and ICD10 code, get the date of the first diagnose

def main():
    verbosity=logging.INFO

    parser=argparse.ArgumentParser(description="Given patient ID and ICD10 code, get the date of the first diagnosis.")
    parser.add_argument('--project','-p',required=True,action='store',help="Project name")
    parser.add_argument('--release','-r',required=True,action='store',help="Project release")
    parser.add_argument('-icd10','--icd10',required=True,action='store',help="ICD10 code")
    parser.add_argument('-id','--id',required=True,action='store',help="Patient ID")
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
    ID=args.id
    icd10=args.icd10

    if args.verbose is not None:
        if args.verbose=="debug":
            verbosity=logging.DEBUG
        elif args.verbose=="warning":
            verbosity=logging.WARNING
        elif args.verbose=="error":
            verbosity=logging.ERROR

    LOGGER=logging.getLogger("first_diagnosis_date")
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
    LOGGER.info("")

#-----------------------------------------------------------------------------------------------------------------------------

    with tarfile.open(infile,"r:*") as tar:
        df_main=pd.read_table(tar.extractfile("hesin.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","ins_index","epistart"])
        df_diag=pd.read_table(tar.extractfile("hesin_diag.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","ins_index","diag_icd10"])
        JT=pd.merge(df_main,df_diag,on=["eid","ins_index"],how="inner")
        JT=JT[(JT["diag_icd10"]==icd10) & (JT["eid"]==ID)]
        if (len(JT)==0):        
            LOGGER.info("Could not find any records for ID=%s, ICD10=%s" %(ID,icd10))
        else:
            JT["epistart_fmt"]=pd.to_datetime(JT["epistart"])
            JT.sort_values(by="epistart_fmt",ascending=True,inplace=True)
            print(JT.head(1).rename(columns={"eid":"ID","epistart":"Date","diag_icd10":"ICD10"}).to_csv(sep="\t",index=False,columns=["ID","Date","ICD10"]),end='')
        
if __name__=="__main__":
    main()
