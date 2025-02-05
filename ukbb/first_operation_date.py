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

# opdate -> epistart -> admidate -> epiend
def calc_operation_date(row):
    if row["opdate"]=="NA":
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
    else:
        return row["opdate"]

def main():
    verbosity=logging.INFO

    parser=argparse.ArgumentParser(description="Given patient IDs and OPCS4 code(s), get the date of the first operation",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    rgroup=parser.add_argument_group("required arguments")
    rgroup.add_argument('--project','-p',required=True,action='store',help="Project name")
    rgroup.add_argument('--release','-r',required=True,action='store',help="Project HESIN release")
    group=rgroup.add_mutually_exclusive_group(required=True)
    group.add_argument('--opcs4','-opcs4',required=False,action='store',help="OPCS4 code or a file/list with OPCS4 codes")
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
    opcs4=utils.readIDList(args.opcs4)

    if args.verbose is not None:
        if args.verbose=="debug":
            verbosity=logging.DEBUG
        elif args.verbose=="warning":
            verbosity=logging.WARNING
        elif args.verbose=="error":
            verbosity=logging.ERROR

    LOGGER=logging.getLogger("first_operation_date")
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
    LOGGER.debug("total OPCS4 codes: %d" %(len(opcs4)))

#-----------------------------------------------------------------------------------------------------------------------------

    with tarfile.open(infile,"r:*") as tar:
        df_main=pd.read_table(tar.extractfile("hesin.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","ins_index","epistart","epiend","admidate"])
        df_oper=pd.read_table(tar.extractfile("hesin_oper.txt"),sep="\t",header=0,dtype=str,quotechar='"',quoting=csv.QUOTE_NONE,keep_default_na=False,usecols=["eid","ins_index","opdate","oper4"])
        JT=pd.merge(df_main,df_oper,on=["eid","ins_index"],how="inner")
        if id_list is None or len(id_list)==0:
            JT=JT[JT["oper4"].isin(opcs4)]
        else:
            JT=JT[(JT["oper4"].isin(opcs4)) & (JT["eid"].isin(id_list))]
                
        # print(JT.to_csv(sep="\t",index=False),end='',file=sys.stderr)
        if (len(JT)==0):
            LOGGER.info("could not find any records for provided OPCS4 codes / sample IDs")
        else:
            JT["operation_date"]=JT.apply(calc_operation_date,axis=1)
            JT["operation_date_fmt"]=pd.to_datetime(JT["operation_date"],dayfirst=True,errors="coerce")
            # only interested in the earliest date, NaT values are ignored
            # if for a sample there are only NaT values then this sample will not be in idx, so not reported
            idx=JT.groupby(["eid"])["operation_date_fmt"].transform(min)==JT["operation_date_fmt"]
            # sometimes for the same ID and ICD there are multiple entries with the same date, so we need drop_duplicates
            if args.output:
                JT[idx][["eid","operation_date"]].drop_duplicates().rename(columns={"eid":"ID","operation_date":"date"}).to_csv(args.output,sep="\t",index=False)
            else:
                print(JT[idx][["eid","operation_date"]].drop_duplicates().rename(columns={"eid":"ID","operation_date":"date"}).to_csv(sep="\t",index=False),end='')
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
