#!/usr/bin/env bash

scriptname=$0
args=("$@")

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
collapsescript="${scriptdir}/collapseFields.pl"
collectstats3="${scriptdir}/collectStats.pl"

function usage () {
    echo ""
    echo "Script for selecting data from a UKBB release"
    echo ""
    echo "Usage: ukbb_select.sh -p | --project <project name>"
    echo "                      -r | --release <release>"
    echo "                      -l | --list-fields <output available fields and exit>"
    echo "                      --majority <comma-separated list of fields>"
    echo "                      --mean <comma-separated list of fields>"
    echo "                      --min-missing <comma-separated list of fields>"
    echo "                      --cc <field,value>"
    echo "                      -o | --output <optional: output prefix; if not specified, output goes to STDOUT>"
    echo "                      -c | --config <optional: config file; default: config.txt in script directory>"
    echo "                      -n | --use-names <optional: use trait names in output header; default: false>"
    echo "                      -h | --help"
    echo ""
    echo "--majority: for each ID output most frequent value across instances/array indexes; useful for categorical variables"
    echo "--mean: for each ID output mean value across instances/array indexes; useful for integer/continuous variables"
    echo "--min-missing: select instance/array index with the least number of NAs"
    echo "--cc: case/control: for each ID output \"1\" if a specified value occurs among field's instances/array indexes, otherwise \"0\""
    echo ""
    exit 0
}

if [[ $# -eq 0 ]];then
    usage
fi

OPTS=$(getopt -o hnlc:p:r:o: -l help,use-names,list-fields,project:,release:,config:,mean:,majority:,min-missing:,cc:,output: -n 'ukbb_select' -- "$@")

if [ $? != 0 ] ; then echo "ERROR: failed parsing options" >&2 ; usage ; exit 1 ; fi

eval set -- "$OPTS"

declare -a ccargs
declare -a majorityargs
declare -a meanargs
declare -a minnaargs

declare -A ccfields # field --> val1,val2,val3,...
declare -A majorityfields # field --> 1
declare -A meanfields # field --> 1
declare -A minnafields # field --> 1

declare -A available_projects

usenames="NO"
config=""
project=""
release=""
outfile=""
use_gzip="NO"
list_fields="NO"
while true; do
  case "$1" in
    -h|--help ) usage; shift ;;
    -n|--use-names ) usenames="YES"; shift ;;
    -l|--list-fields ) list_fields="YES"; shift ;;
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

if [[ ! -z "$outfile" ]];then
    logfile="$outfile".log
    outfile="$outfile".txt.gz
    exitIfExists "$outfile" "ERROR: output file $outfile already exists"
    : > "$logfile"
    sfx="|tee -a $logfile"
else
    sfx="1>&2"
fi

eval date "+%d-%b-%Y:%H-%M-%S" "$sfx"
eval echo "" "$sfx"
eval echo "Current dir: ${PWD}" "$sfx"
eval echo "Command line: $scriptname ${args[@]}" "$sfx"
eval echo "" "$sfx"

if [[ "$outfile" =~ \.gz$ ]];then
    use_gzip="YES"    
fi

if [[ -z "$config" ]];then
    config="${scriptdir}"/config.txt
    eval echo "INFO: no config file specified, using $config" "$sfx"
fi
exitIfNotFile "$config" "ERROR: config $config does not exist"

readValue "$config" DATA_DICT dfile
exitIfNotFile "$dfile" "ERROR: dictionary $dfile does not exist"
eval echo "INFO: using data dictionary: $dfile" "$sfx"

readAArray "$config" PROJECTS available_projects
exitIfEmpty "${available_projects[$project]}" "ERROR: project $project is not defined in $config"
readValue "$config" DATA_PATH data_path
exitIfNotDir "$data_path" "ERROR: $data_path is not a directory"

data_path=${data_path%/}
infile="$data_path"/"${available_projects[$project]}"/releases/phenotypes_r"$release".txt.gz
eval echo "INFO: using input file $infile" "$sfx"
exitIfNotFile "$infile" "ERROR: input file $infile does not exist"

fieldID_cn=$(getColNum "$dfile" "FieldID" "cat")
field_cn=$(getColNum "$dfile" "Field" "cat")
valtype_cn=$(getColNum "$dfile" "ValueType" "cat")
ID_cn=$(getColNum "$infile" "f.eid" "zcat")

eval echo "INFO: field ID column in data dictionary: $fieldID_cn" "$sfx"
eval echo "INFO: field description column in data dictionary: $field_cn" "$sfx"
eval echo "INFO: value type column in data dictionary: $valtype_cn" "$sfx"
eval echo "INFO: ID column in $infile: $ID_cn" "$sfx"

# list fields and exit
if [[ "$list_fields" == "YES" ]];then
    declare -a colnames
    declare -A h1
    getColNames "$infile" "zcat" colnames
    for c in "${colnames[@]}";do
	n=$(tail -n +2 "$dfile"| awk -v n1=$fieldID_cn -v n2=$field_cn -v n3=$c 'BEGIN{FS=OFS="\t";ret="NA";}$n1==n3{ret=$n2;exit;}END{print ret;}')
	h1["$c"]="$n"
    done
    if [[ ! -z "$outfile" ]];then
	if [[ "$use_gzip" == "NO" ]];then
	    for c in "${!h1[@]}";do
		echo -e "$c\t${h1[$c]}"
	    done | sort -k1,1n > "$outfile"
	else
	    for c in "${!h1[@]}";do
		echo -e "$c\t${h1[$c]}"
	    done | sort -k1,1n | gzip - -c > "$outfile"
	fi
    else
	for c in "${!h1[@]}";do
	    echo -e "$c\t${h1[$c]}"
	done | sort -k1,1n 
    fi
    
    eval date "+%d-%b-%Y:%H-%M-%S" "$sfx"
    exit 0
fi

if [[ "$usenames" == "YES" ]];then
    eval echo "INFO: using human readable field names in output" "$sfx"
fi

# check if input fields are integers
# and
# filling *fields associative arrays
eval echo "" "$sfx"
eval echo "INFO: checking if input fields are integers" "$sfx"
for name in mean minna majority;do
    declare -n Z=${name}args
    declare -n Y=${name}fields
    for x in "${Z[@]}";do
	while read z;do
	    if [[ $z =~ ^[0-9]+$ ]];then
		Y[$z]=1
	    else
		eval echo 'WARN: field \"'"$z"'\" is not an integer\; dropping' "$sfx"
	    fi
	done < <(echo $x|tr ',' '\n')
    done
done

declare -A temp_ar
for x in "${ccargs[@]}";do
    read a b < <(echo $x|tr ',' ' ')
    if [[ $a =~ ^[0-9]+$ ]];then
	v="${ccfields[$a]}"
	if [[ -z "$v" ]];then
	    ccfields[$a]=$b
	else
	    temp_ar=()
	    for x in $(echo "$v"|tr ',' ' '); do temp_ar[$x]=1; done
	    # only add if not already present
	    z="${temp_ar[$b]}"
	    if [[ -z "$z" ]];then
		ccfields[$a]=$v","$b
	    fi
	fi
    else
	eval echo 'WARN: specified field \"'"$a"'\" is not an integer\; dropping' "$sfx"
    fi		
done
eval echo "" "$sfx"

# check if input fields are present
declare -a to_delete
declare -A ar
eval echo "INFO: checking if input fields are present in input" "$sfx"
for name in ccfields meanfields minnafields majorityfields;do
    declare -n Z=$name
    for f in "${!Z[@]}";do
	getCols "$infile" "zcat" "$f" ar
	if [[ "${#ar[@]}" -eq 0 ]];then
	    eval echo "INFO: checking $f: not in input $sfx"
	    to_delete+=($f)
	else
	    eval echo "INFO: checking $f: OK $sfx"
	fi
    done
    for f in "${to_delete[@]}";do
	unset "Z[$f]"
    done
    to_delete=()
done
eval echo "" "$sfx"

# checking field types
declare -A fdesc
declare -A ftype

eval echo "INFO: getting field value type information" "$sfx"
while IFS=$'\t' read fid f t;do
    fdesc["$fid"]="$f"
    ftype["$fid"]="$t"
done < <(tail -n +2 "$dfile"| awk -v n1=$fieldID_cn -v n2=$field_cn -v n3=$valtype_cn 'BEGIN{FS=OFS="\t";}{print $n1,$n2,$n3;}')

to_delete=()
for f in "${!majorityfields[@]}";do
    t="${ftype[$f]}"
    if [[ -z "$t" ]];then
	eval echo 'WARN: no value type information found for '"$f"'\; dropping' "$sfx"
	to_delete+=($f)
    else
	if [[ "$t" != "Categorical single" ]];then
	    eval echo 'WARN: '"$f"' has value type \"'"$t"'\"\; for the majority rule it needs to be \"Categorical single\"\; dropping' "$sfx"
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
	eval echo 'WARN: no value type information found for '"$f"'\; dropping' "$sfx"
	to_delete+=($f)
    else
	if [[ ! "$t" =~ ^Categorical ]];then
	    eval echo 'WARN: '"$f"' has value type \"'"$t"'\"\; for case-control it needs to be \"Categorical single\" or \"Categorical multiple\"\; dropping' "$sfx"
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
	eval echo 'WARN: no value type information found for '"$f"'\; dropping' "$sfx"
	to_delete+=($f)
    else
	if [[ "$t" != "Continuous" && "$t" != "Integer" ]];then
	    eval echo 'WARN: '"$f"' has value type \"'"$t"'\"\; for the mean rule it needs to be \"Continuous\" or \"Integer\"\; dropping' "$sfx"
	    to_delete+=($f)
	fi
    fi
done
for f in "${to_delete[@]}";do
    unset "meanfields[$f]"
done
to_delete=()
eval echo "" "$sfx"

# ----------------------------------------------------------------------
# reporting

if [[ "${#ccfields[@]}" -eq 0 && "${#meanfields[@]}" -eq 0 && "${#majorityfields[@]}" -eq 0 && "${#minnafields[@]}" -eq 0 ]];then
    eval echo 'INFO: no fields left to process; exit' "$sfx"
    eval echo "" "$sfx"
    eval date "+%d-%b-%Y:%H-%M-%S" "$sfx"
    exit 0
fi

if [[ "${#ccfields[@]}" -ne 0 ]];then
    eval echo 'INFO: selecting case/control fields:' "$sfx"
    for f in "${!ccfields[@]}";do
	v="${ccfields[$f]}"
	eval echo "$f : values $v" "$sfx"
    done
    eval echo "" "$sfx"
fi

if [[ "${#meanfields[@]}" -ne 0 ]];then
    eval echo 'INFO: selecting mean value fields:' "$sfx"
    for f in "${!meanfields[@]}";do
	eval echo "$f" "$sfx"
    done
    eval echo "" "$sfx"
fi

if [[ "${#majorityfields[@]}" -ne 0 ]];then
    eval echo 'INFO: selecting majority rule fields:' "$sfx"
    for f in "${!majorityfields[@]}";do
	eval echo "$f" "$sfx"
    done
    eval echo "" "$sfx"
fi

if [[ "${#minnafields[@]}" -ne 0 ]];then
    eval echo 'INFO: selecting min NA fields:' "$sfx"
    for f in "${!minnafields[@]}";do
	eval echo "$f" "$sfx"
    done
    eval echo "" "$sfx"
fi

# ----------------------------------------------------------------------

# selecting
outdir=$(dirname "$outfile")
tmpdir=$(mktemp -d -p "$outdir" tempdir_select_XXXXXXXX)
if [[ $? -ne 0 ]];then
    eval echo "ERROR: could not create temporary directory in $outdir " "$sfx"
    exit 1
fi

eval echo "INFO: selecting fields" "$sfx"
for f in "${!majorityfields[@]}";do
    getCols "$infile" "zcat" "$f" ar
    str=$(join_by , "${!ar[@]}")
    colname="$f"
    if [[ "$usenames" == "YES" ]];then
	if [[ ! -z "${fdesc[$f]}" ]];then
	    colname="${fdesc[$f]}"
	fi
    fi
    colname="$colname (majority)"
    echo -e "f.eid\t$colname" > "$tmpdir"/majority_"$f"
    paste <(zcat "$infile"|cut -f "$ID_cn"|tail -n +3) <(zcat "$infile"|cut -f "$str"|tail -n +3)|"$collapsescript" majority >> "$tmpdir"/majority_"$f"
done

for f in "${!meanfields[@]}";do
    getCols "$infile" "zcat" "$f" ar
    str=$(join_by , "${!ar[@]}")
    d="${fdesc[$f]}"
    colname="$f"
    if [[ "$usenames" == "YES" ]];then
	if [[ ! -z "$d" ]];then
	    colname="$d"
	fi
    fi
    colname="$colname (mean)"
    echo -e "f.eid\t$colname" > "$tmpdir"/mean_"$f"
    paste <(zcat "$infile"|cut -f "$ID_cn"|tail -n +3) <(zcat "$infile"|cut -f "$str"|tail -n +3)|"$collapsescript" mean >> "$tmpdir"/mean_"$f"
done

for f in "${!ccfields[@]}";do
    getCols "$infile" "zcat" "$f" ar
    str=$(join_by , "${!ar[@]}")
    d="${fdesc[$f]}"
    val="${ccfields[$f]}"
    for v in $(echo "${val}"|tr ',' ' '); do
	colname="$f:${v}"
	if [[ "$usenames" == "YES" ]];then
	    if [[ ! -z "$d" ]];then
		colname="$d:${v}"
	    fi
	fi
	colname="$colname (case/control)"
	echo -e "f.eid\t$colname" > "$tmpdir"/cc_"${f}"_"${v}"
	paste <(zcat "$infile"|cut -f "$ID_cn"|tail -n +3) <(zcat "$infile"|cut -f "$str"|tail -n +3)|"$collapsescript" cc ${v} >> "$tmpdir"/cc_"${f}"_"${v}"
    done
done

for f in "${!minnafields[@]}";do
    getCols "$infile" "zcat" "$f" ar
    str=$(join_by , "${!ar[@]}")
    f1=""
    if [[ "${#ar[@]}" -eq 1 ]];then
	f1="${ar[$str]}"
	eval echo "INFO: there is only one column in input for field $f: $f1" "$sfx"
    else
	f1=$(zcat "$infile"|cut -f "$str"|"$collectstats3"|sort -k2,2n|head -n 1|cut -f 1)	
	eval echo "INFO: selecting column with least NAs for field $f: $f1" "$sfx"
    fi
    colname="${f1}"
    d="${fdesc[$f]}"
    if [[ "$usenames" == "YES" ]];then
	if [[ ! -z "$d" ]];then
	    colname="$d"
	fi
    fi
    colname="$colname (min NA)"
    echo -e "f.eid\t$colname" > "$tmpdir"/minna_"$f"
    n=$(getColNum "$infile" "${f1}" "zcat")
    paste <(zcat "$infile"|cut -f "$ID_cn"|tail -n +3) <(zcat "$infile"|cut -f "$n"|tail -n +3) >> "$tmpdir"/minna_"$f"
done
eval echo "" "$sfx"

# merging
eval echo "INFO: merging" "$sfx"
i=0
for f in $(find "$tmpdir" -maxdepth 1 -type f);do
    if [[ "$i" -eq 0 ]];then
	cmd="paste $f"
    else
	cmd="$cmd <(cut -f 2 $f)"
    fi
    i=$((i+1))
done

if [[ ! -z "$outfile" ]];then
    if [[ "$use_gzip" == "NO" ]];then
	eval "$cmd" > "$outfile"
    else
	eval "$cmd" | gzip - -c > "$outfile"
    fi
else
    eval "$cmd"
fi

eval echo "" "$sfx"
eval date "+%d-%b-%Y:%H-%M-%S" "$sfx"

rm -rf "$tmpdir"
exit 0



	 
