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

# given patient ID(s) and ICD10 code, get the date of the first diagnosis

# epistart -> admidate -> epiend
def calc_diagnosis_date(row):
    if row["epistart"]=="NA":
        if row["admidate"]=="NA":
            if row["epiend"]=="NA":
                return "NA"
            else:
                return row["epiend"]
        else:
            return row["admidate"]
    else:
        return row["epistart"]

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

    parser=argparse.ArgumentParser(description="Given patient ID and ICD10 code, get the date of the first diagnosis.")
    parser.add_argument('--project','-p',required=True,action='store',help="Project name")
    parser.add_argument('--release','-r',required=True,action='store',help="Project release")
    parser.add_argument('--icd10','-icd10',required=True,action='store',help="ICD10 code")
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

    # if len(id_list)==0:
    #     LOGGER.error("No patient IDs provided")
    #     sys.exit(1)
    # else:
    #     LOGGER.info(str(len(id_list))+" patient ID(s) provided")

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

#-----------------------------------------------------------------------------------------------------------------------------

    with tarfile.open(infile,"r:*") as tar:
        df_main=pd.read_table(tar.extractfile("hesin.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","ins_index","epistart","epiend","admidate"])
        df_diag=pd.read_table(tar.extractfile("hesin_diag.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","ins_index","diag_icd10"])
        JT=pd.merge(df_main,df_diag,on=["eid","ins_index"],how="inner")
        if id_list is None or len(id_list)==0:
            JT=JT[JT["diag_icd10"]==icd10]
        else:
            JT=JT[(JT["diag_icd10"]==icd10) & (JT["eid"].isin(id_list))]
        # print(JT.to_csv(sep="\t",index=False),end='',file=sys.stderr)
        if (len(JT)==0):        
            LOGGER.info("Could not find any records for ICD10=%s" %(icd10))
        else:
            JT["diagnosis_date"]=JT.apply(calc_diagnosis_date,axis=1)
            JT["diagnosis_date_fmt"]=pd.to_datetime(JT["diagnosis_date"],dayfirst=True,errors="coerce")
            # only interested in the earliest date, NaT values are ignored
            # if for a sample there are only NaT values then this sample will not be in idx, so not reported
            idx=JT.groupby(["eid"])["diagnosis_date_fmt"].transform(min)==JT["diagnosis_date_fmt"]
            # sometimes for the same ID and ICD10 there are multiple entries with the same date, so we need drop_duplicates
            print(JT[idx][["eid","diagnosis_date","diag_icd10"]].drop_duplicates().rename(columns={"eid":"ID","diagnosis_date":"Date","diag_icd10":"ICD10"}).to_csv(sep="\t",index=False),end='')
        
if __name__=="__main__":
    main()
