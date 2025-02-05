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

# given patient ID(s) and ICD code, get the date of the first diagnosis

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

def main():
    verbosity=logging.INFO

    parser=argparse.ArgumentParser(description="Given patient IDs and ICD code(s), get the date of the first diagnosis",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    rgroup=parser.add_argument_group("required arguments")
    rgroup.add_argument('--project','-p',required=True,action='store',help="Project name")
    rgroup.add_argument('--release','-r',required=True,action='store',help="Project HESIN release")
    group=rgroup.add_mutually_exclusive_group(required=True)
    group.add_argument('--icd10','-icd10',required=False,action='store',help="ICD10 code or a file/list with ICD10 codes")
    group.add_argument('--icd9','-icd9',required=False,action='store',help="ICD9 code or a file/list with ICD9 codes")
    # opgroup=parser.add_argument_group("optional arguments")
    parser.add_argument('--id','-id',required=False,action='store',help="Patient ID or a file/list with patient IDs")
    parser.add_argument('--output','-o',required=False,action='store',help="Output file")
    parser.add_argument('--config','-c',required=False,action='store',help="Config file")
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
    LOGGER.debug("config file: %s" % config)
    CONFIG=utils.readConfig(config)
    if CONFIG is None:
        sys.exit(1)
    infile=utils.getProjectFileName(CONFIG,project,release,"HESIN")
    if infile is None:
        sys.exit(1)

    LOGGER.debug("input file: %s" % infile)
    if id_list is None or len(id_list)==0:
        LOGGER.info("using all samples in input file")
    LOGGER.info("")

    icd=icd9
    icd_col="diag_icd9"
    if icd9 is None:
        icd=icd10
        icd_col="diag_icd10"

    LOGGER.debug("total ICD9 codes: %d" %(len(icd9) if icd9 else 0))
    LOGGER.debug("total ICD10 codes: %d" %(len(icd10) if icd10 else 0))

#-----------------------------------------------------------------------------------------------------------------------------

    with tarfile.open(infile,"r:*") as tar:
        df_main=pd.read_table(tar.extractfile("hesin.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","ins_index","epistart","epiend","admidate"])
        df_diag=pd.read_table(tar.extractfile("hesin_diag.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","ins_index",icd_col])
        JT=pd.merge(df_main,df_diag,on=["eid","ins_index"],how="inner")
        if id_list is None or len(id_list)==0:
            JT=JT[JT[icd_col].isin(icd)]
        else:
            JT=JT[(JT[icd_col].isin(icd)) & (JT["eid"].isin(id_list))]
                
        # print(JT.to_csv(sep="\t",index=False),end='',file=sys.stderr)
        if (len(JT)==0):
            LOGGER.info("could not find any records for provided ICD codes / sample IDs")
        else:
            JT["diagnosis_date"]=JT.apply(calc_diagnosis_date,axis=1)
            JT["diagnosis_date_fmt"]=pd.to_datetime(JT["diagnosis_date"],dayfirst=True,errors="coerce")
            # only interested in the earliest date, NaT values are ignored
            # if for a sample there are only NaT values then this sample will not be in idx, so not reported
            idx=JT.groupby(["eid"])["diagnosis_date_fmt"].transform(min)==JT["diagnosis_date_fmt"]
            # sometimes for the same ID and ICD there are multiple entries with the same date, so we need drop_duplicates
            if args.output:
                JT[idx][["eid","diagnosis_date"]].drop_duplicates().rename(columns={"eid":"ID","diagnosis_date":"date"}).to_csv(args.output,sep="\t",index=False)
            else:
                print(JT[idx][["eid","diagnosis_date"]].drop_duplicates().rename(columns={"eid":"ID","diagnosis_date":"date"}).to_csv(sep="\t",index=False),end='')
            if not(id_list is None) and len(id_list)!=0:
                # IDs in input list but not in the previous output
                not_found_ids=not_found_ids.union(set(id_list)-set(JT[idx]["eid"]))
                if len(not_found_ids)!=0:
                    if args.output:
                        pd.DataFrame({"A":list(not_found_ids),"B":"NA"}).to_csv(args.output,header=False,index=False,sep="\t")
                    else:
                        print(pd.DataFrame({"A":list(not_found_ids),"B":"NA"}).to_csv(header=False,index=False,sep="\t"),end='')
                        
if __name__=="__main__":
    main()
