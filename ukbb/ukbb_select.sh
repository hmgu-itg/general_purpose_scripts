#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
collapsescript="${scriptdir}/collapseFields.pl"
collectstats3="${scriptdir}/collectStats3.pl"

OPTS=$(getopt -o hnc:p:r:o: -l help,names,project:,release:,config:,mean:,majority:,min-missing:,cc:,output: -n 'ukbb_select' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

declare -a ccargs
declare -a majorityargs
declare -a meanargs
declare -a minnaargs

declare -A ccfields
declare -A majorityfields
declare -A meanfields
declare -A minnafields

declare -A available_projects

usenames="NO"
config=""
project=""
release=""
outfile=""
use_gzip="NO"
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
    -n|--names ) usenames="YES"; shift ;;
    -p|--project ) project=$2; shift 2 ;;
    -r|--release ) release=$2; shift 2 ;;
    -c|--config ) config=$2; shift 2 ;;
    --majority ) majorityargs+=($2); shift 2 ;;
    --mean ) meanargs+=($2); shift 2 ;;
    --min-missing ) minnaargs+=($2); shift 2 ;;
    --cc ) ccargs+=($2); shift 2 ;;
    -o|--output ) outfile=$2; shift 2 ;;
    --) shift ; break ;;
    * ) echo "Unexpected option: $1" ; usage ;;
  esac
done

exitIfEmpty "$project" "ERROR: project not specified"
exitIfEmpty "$release" "ERROR: release not specified"

logfile="$outfile".log
outfile="$outfile".txt.gz
exitIfEmpty "$outfile" "ERROR: output file not specified"
exitIfExists "$outfile" "ERROR: output file $outfile already exists"

: > "$logfile"
date "+%F %H-%M-%S"|tee -a "$logfile"
echo "Current dir: ${PWD}"|tee -a "$logfile"
echo "Command line: $scriptname ${args[@]}"|tee -a "$logfile"
echo ""|tee -a "$logfile"

if [[ "$outfile" =~ \.gz$ ]];then
    use_gzip="YES"    
fi

if [[ -z "$config" ]];then
    config="${scriptdir}"/config.txt
    echo "INFO: no config file specified, using $config"|tee -a "$logfile"
fi
exitIfNotFile "$config" "ERROR: config $config does not exist"

readValue "$config" DATA_DICT dfile
exitIfNotFile "$dfile" "ERROR: dictionary $dfile does not exist"

readAArray "$config" PROJECTS available_projects
exitIfEmpty "${available_projects[$project]}" "ERROR: project $project is not defined in $config"
readValue "$config" DATA_PATH data_path
exitIfNotDir "$data_path" "ERROR: $data_path is not a directory"

data_path=${data_path%/}
infile="$data_path"/"${available_projects[$project]}"/releases/phenotypes_r"$release".txt.gz
echo "INFO: using input file $infile"|tee -a "$logfile"
exitIfNotFile "$infile" "ERROR: input file $infile does not exist"

# declare -p ccargs
# declare -p majorityargs
# declare -p meanargs
# declare -p minnaargs

echo "INFO: using human readable field names in output: $usenames"|tee -a "$logfile"

for name in mean minna majority;do
    declare -n Z=${name}args
    declare -n Y=${name}fields
    for x in "${Z[@]}";do
	while read z;do
	    if [[ $z =~ ^[0-9]+$ ]];then
		Y[$z]=1
	    else
		echo "WARN: specified field $z is not an integer; dropping"|tee -a "$logfile"
	    fi
	done < <(echo $x|tr ',' '\n')
    done
done

for x in "${ccargs[@]}";do
    read a b < <(echo $x|tr ',' ' ')
    if [[ $a =~ ^[0-9]+$ ]];then
	ccfields[$a]=$b
    else
	echo "WARN: specified field $a is not an integer; dropping"|tee -a "$logfile"
    fi		
done

# declare -p ccfields
# declare -p majorityfields
# declare -p meanfields
# declare -p minnafields

fieldID_cn=$(getColNum "$dfile" "FieldID" "cat")
field_cn=$(getColNum "$dfile" "Field" "cat")
valtype_cn=$(getColNum "$dfile" "ValueType" "cat")
ID_cn=$(getColNum "$infile" "f.eid" "zcat")

declare -a to_delete
declare -A ar
for name in ccfields meanfields minnafields majorityfields;do
    declare -n Z=$name
    for f in "${!Z[@]}";do
	echo -n "INFO: checking $f: "|tee -a "$logfile"
	getCols "$infile" "zcat" "$f" ar
	if [[ "${#ar[@]}" -eq 0 ]];then
	    echo "not in input"|tee -a "$logfile"
	    to_delete+=($f)
	else
	    echo "OK"|tee -a "$logfile"
	fi
    done
    for f in "${to_delete[@]}";do
	unset "Z[$f]"
    done
    to_delete=()
done

# declare -p ccfields
# declare -p majorityfields
# declare -p meanfields
# declare -p minnafields

# checking field types
declare -A fdesc
declare -A ftype

echo "INFO: getting value type information"|tee -a "$logfile"
while IFS=$'\t' read fid f t;do
    fdesc["$fid"]="$f"
    ftype["$fid"]="$t"
done < <(tail -n +2 "$dfile"| awk -v n1=$fieldID_cn -v n2=$field_cn -v n3=$valtype_cn 'BEGIN{FS=OFS="\t";}{print $n1,$n2,$n3;}')

# fid="21"
# echo "$fid" "${fdesc[$fid]}"
# echo "$fid" "${fdesc[21]}"
# x="${fdesc[$fid]}"
# echo "$fid" "$x"
# if [[ ! -z "${fdesc[$fid]}" ]];then
#     echo "NOT EMPTY: ${fdesc[$fid]}"
# else
#     echo "EMPTY: ${fdesc[$fid]}"
# fi
# exit 1

#declare -p fdesc

to_delete=()
for f in "${!majorityfields[@]}";do
    t="${ftype[$f]}"
    if [[ -z "$t" ]];then
	echo "WARN: no value type information found for $f; dropping"|tee -a "$logfile"
	to_delete+=($f)
    else
	if [[ "$t" != "Categorical single" ]];then
	    echo "WARN: $f has value type \"$t\"; for majority rule it needs to be \"Categorical single\"; dropping"|tee -a "$logfile"
	    to_delete+=($f)
	fi
    fi
done
for f in "${to_delete[@]}";do
    unset "majorityfields[$f]"
done
to_delete=()

for f in "${!ccfields[@]}";do
    t="${ftype[$f]}"
    if [[ -z "$t" ]];then
	echo "WARN: no value type information found for $f; dropping"|tee -a "$logfile"
	to_delete+=($f)
    else
	if [[ ! "$t" =~ ^Categorical ]];then
	    echo "WARN: $f has value type \"$t\"; for case-control it needs to be \"Categorical single\" or \"Categorical multiple\"; dropping"|tee -a "$logfile"
	    to_delete+=($f)
	fi
    fi
done
for f in "${to_delete[@]}";do
    unset "ccfields[$f]"
done
to_delete=()

for f in "${!meanfields[@]}";do
    t="${ftype[$f]}"
    if [[ -z "$t" ]];then
	echo "WARN: no value type information found for $f; dropping"|tee -a "$logfile"
	to_delete+=($f)
    else
	if [[ "$t" != "Continuous" || "$t" != "Integer" ]];then
	    echo "WARN: $f has value type \"$t\"; for mean rule it needs to be \"Continuous\" or \"Integer\"; dropping"|tee -a "$logfile"
	    to_delete+=($f)
	fi
    fi
done
for f in "${to_delete[@]}";do
    unset "meanfields[$f]"
done
to_delete=()

# declare -p ccfields
# declare -p majorityfields
# declare -p meanfields
# declare -p minnafields

outdir=$(dirname "$outfile")
tmpdir=$(mktemp -d -p "$outdir" tempdir_select_XXXXXXXX)
if [[ $? -ne 0 ]];then
    echo "ERROR: could not create temporary directory in $outdir "|tee -a "$logfile"
    exit 1
fi

for f in "${!majorityfields[@]}";do
    getCols "$infile" "zcat" "$f" ar
    #declare -p ar
    str=$(join_by , "${!ar[@]}")
    #echo "str=$str"
    
    # if [[ ! -z "${fdesc[$f]}" ]];then
    # 	echo "NOT EMPTY: ${fdesc[$f]}"
    # else
    # 	echo "EMPTY: ${fdesc[$f]}"
    # fi
	
    # d="${fdesc[$f]}"
    colname="$f"
    # if [[ "$usenames" == "YES" ]];then
    # 	if [[ ! -z "$d" ]];then
    # 	    colname="$d"
    # 	fi
    # fi
    if [[ "$usenames" == "YES" ]];then
	if [[ ! -z "${fdesc[$f]}" ]];then
	    colname="${fdesc[$f]}"
	fi
    fi
    #echo "colname=$colname"
    echo -e "f.eid\t$colname" > "$tmpdir"/majority_"$f"
    paste <(zcat "$infile"|cut -f "$ID_cn"|tail -n +3) <(zcat "$infile"|cut -f "$str"|tail -n +3)|"$collapsescript" majority >> "$tmpdir"/majority_"$f"
done

for f in "${!meanfields[@]}";do
    getCols "$infile" "zcat" "$f" ar
    #declare -p ar
    str=$(join_by , "${!ar[@]}")
    #echo "str=$str"
    d="${fdesc[$f]}"
    colname="$f"
    if [[ "$usenames" == "YES" ]];then
	if [[ ! -z "$d" ]];then
	    colname="$d"
	fi
    fi
    echo -e "f.eid\t$colname" > "$tmpdir"/mean_"$f"
    paste <(zcat "$infile"|cut -f "$ID_cn"|tail -n +3) <(zcat "$infile"|cut -f "$str"|tail -n +3)|"$collapsescript" mean >> "$tmpdir"/mean_"$f"
done

for f in "${!ccfields[@]}";do
    getCols "$infile" "zcat" "$f" ar
    val="${ccfields[$f]}"
    #declare -p ar
    str=$(join_by , "${!ar[@]}")
    #echo "str=$str"
    colname="$f : ${ccfields[$f]}"
    d="${fdesc[$f]}"
    if [[ "$usenames" == "YES" ]];then
	if [[ ! -z "$d" ]];then
	    colname="$d : ${ccfields[$f]}"
	fi
    fi    
    echo -e "f.eid\t$colname" > "$tmpdir"/cc_"$f"
    paste <(zcat "$infile"|cut -f "$ID_cn"|tail -n +3) <(zcat "$infile"|cut -f "$str"|tail -n +3)|"$collapsescript" cc $val >> "$tmpdir"/cc_"$f"
done

for f in "${!minnafields[@]}";do
    getCols "$infile" "zcat" "$f" ar
    #declare -p ar
    str=$(join_by , "${!ar[@]}")
    #echo "str=$str"
    f1=""
    if [[ "${#ar[@]}" -eq 1 ]];then
	f1="${ar[$str]}"
	echo "INFO: there is only one column in input for field $f: $f1"|tee -a "$logfile"
    else
	f1=$(zcat "$infile"|cut -f "$str"|"$collectstats3"|sort -k2,2n|head -n 1|cut -f 1)	
	echo "INFO: selecting column with least NAs for field $f: $f1"|tee -a "$logfile"
    fi
    #echo "f1=${f1}"
    colname="${f1}"
    d="${fdesc[$f]}"
    if [[ "$usenames" == "YES" ]];then
	if [[ ! -z "$d" ]];then
	    colname="$d"
	fi
    fi    
    echo -e "f.eid\t$colname" > "$tmpdir"/minna_"$f"
    n=$(getColNum "$infile" "${f1}" "zcat")
    paste <(zcat "$infile"|cut -f "$ID_cn"|tail -n +3) <(zcat "$infile"|cut -f "$n"|tail -n +3) >> "$tmpdir"/minna_"$f"
done

# merging
i=0
for f in $(find "$tmpdir" -maxdepth 1 -type f);do
    if [[ "$i" -eq 0 ]];then
	cmd="paste $f"
    else
	cmd="$cmd <(cut -f 2 $f)"
    fi
    i=$((i+1))
done

if [[ "$use_gzip" == "NO" ]];then
    eval "$cmd" > "$outfile"
else
    eval "$cmd" | gzip - -c > "$outfile"
fi

date "+%F %H-%M-%S"|tee -a "$logfile"

rm -rf "$tmpdir"
exit 0



	 
