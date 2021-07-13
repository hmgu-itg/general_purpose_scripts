#!/usr/bin/env bash

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"

OPTS=$(getopt -o hn -l mean:,median:,majority:,min:,max:,min-missing:,cc: -n 'ukbb_select' -- "$@")

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

declare -a ccargs
declare -a majorityargs
declare -a medianargs
declare -a meanargs
declare -a minargs
declare -a maxargs
declare -a minnaargs
usenames="NO"
meanfields=""
while true; do
  case "$1" in
    -h ) usage ;;
    -n ) usenames="YES" ;;
    --majority ) majorityargs+=$2; shift 2 ;;
    --mean ) meanargs+=$2; shift 2 ;;
    --median ) medianargs+=$2; shift 2 ;;
    --min ) minargs+=$2; shift 2 ;;
    --max ) maxargs+=$2; shift 2 ;;
    --min-missing ) minnaargs+=$2; shift 2 ;;
    --cc ) ccargs+=$2; shift 2 ;;
    --) shift ; break ;;
    * ) echo "Unexpected option: $1" ; usage ;;
  esac
done

declare -p ccargs
declare -p majorityargs
declare -p medianargs
declare -p meanargs
declare -p minargs
declare -p maxargs
declare -p minnaargs
