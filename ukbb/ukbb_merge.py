#!/usr/bin/python3

import argparse
import pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument('-i','--input',required=True,action='store',help="")
parser.add_argument('-u','--update',required=True,action='append',help="")

args=parser.parse_args()

infile=args.input
updates=args.update

print("input: %s updates: %s" %(infile,",".join(updates)))
