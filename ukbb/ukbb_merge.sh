#!/usr/bin/env bash

# ALL FILES ARE TAB SEPARATED

scriptname=$0
args=("$@")

# help folrmatting
bold=$(tput bold)
underlined=$(tput smul)
normal=$(tput sgr0)

scriptdir=$(dirname $(readlink -f $0))
source "${scriptdir}/functions.sh"
collectstats="${scriptdir}/collectStats.pl"
collectstats2="${scriptdir}/collectStats2.pl"

# TODO: check if scripts exist

function usage () {
    echo ""
    echo "Merging/updating script for UKBB data"
    echo ""
    echo "Usage: ukbb_merge.sh -f <ID field name; default: \"f.eid\">"
    echo "                     -i input1.tab"
    echo "                     -i input2.tab"
    echo "           ...                    "
    echo "                     -i inputN.tab"
    echo ""
    echo "                     -u update1.tab"
    echo "                     -u update2.tab"
    echo "           ...                    "
    echo "                     -u updateM.tab"
    echo ""
    echo "                     -o <output dir>"
    echo "                     -r <release; default: \"1\" if merging, incrementing RELEASE in input if updating>"
    echo "                     -k <when updating, also include samples from update file(s) that are ${bold}not${normal} present in input; default: false>"
    echo "                     -x <list of individual IDs to exclude>"
    echo "                     -d <debug mode: do not remove temporary files>"
    echo "                     -s <skip creating missingness summary file>"
    echo ""
    echo "All input/update/output files are tab-separated"
    echo ""
    echo "${underlined}Merge mode${normal}: if no update (-u) files are specified, the script merges all input files."
    echo ""
    echo "${underlined}Update mode${normal}: if at least one update (-u) file is given, the script only uses the first input (-i) file" 
    echo "and ignores the remaining input files. If a column in the input file"
    echo "is present in an update file, its content in the input file will be updated."
    echo ""
    exit 0
}

declare -a input_fnames
declare -a update_fnames
declare -a input_ID_column
declare -a update_ID_column
declare -a input_nrows
declare -a update_nrows
declare -a input_ncols
declare -a update_ncols
declare -a class_array
declare -A column_class
declare -A cats

if [[ $# -eq 0 ]];then
    usage
fi

id_field="f.eid"
datestr=$(date +%F)
out_dir=""
release=""
exclude_list=""
mode="inner"
debug="NO"
skipmiss="NO"
while getopts "hi:u:f:o:r:x:kds" opt; do
    case $opt in
        i)input_fnames+=($OPTARG);;
        u)update_fnames+=($OPTARG);;
        f)id_field=($OPTARG);;
        o)out_dir=($OPTARG);;
        r)release=($OPTARG);;
        x)exclude_list=($OPTARG);;
        k)mode="right";;
        d)debug="YES";;
        s)skipmiss="YES";;
        h)usage;;
        *)usage;;
    esac
done
shift "$((OPTIND-1))"

ret=$(checkInstalledCommands "gawk" "parallel")
if [[ $ret -ne 0 ]];then
    exit 1
fi

if [[ -z "$out_dir" ]];then
    echo "ERROR: no output dir specified" 1>&2
    exit 1
fi

if [[ ! -d "$out_dir" ]];then
    echo "ERROR: output dir $out_dir is not a directory" 1>&2
    exit 1
fi

if [[ ! -w "$out_dir" ]];then
    echo "ERROR: output dir $out_dir is not writable" 1>&2
    exit 1
fi

if [[ ! -z "$exclude_list" ]];then
    if [[ ! -f "$exclude_list" ]];then
	echo "ERROR: exclude list $exclude_list is not a file" 1>&2
	exit 1
    fi
fi

out_dir=${out_dir%/}

n_input="${#input_fnames[@]}"
n_update="${#update_fnames[@]}"

if [[ $n_input -eq 0 ]];then
    echo "ERROR: no input files specified" 1>&2
    exit 1
fi

if [[ $n_input -eq 1 && $n_update -eq 0 && -z "$exclude_list" ]];then
    echo "ERROR: no update files, only one input file specified, no exclude ID list; nothing to do" 1>&2
    exit 1
fi

if [[ $n_input -eq 0 && $n_update -gt 0 ]];then
    echo "ERROR: update file(s) specified but no input files specified" 1>&2
    exit 1
fi

if [[ $n_update -gt 0 && $n_input -gt 1 ]];then
    echo "WARN: $n_update update files and $n_input input files specified; only the first input file (${input_fnames[0]}) will be processed" 1>&2
fi

# get zcat/cat command
for i in $(seq 0 $((n_input-1)));do
    cats["${input_fnames[$i]}"]=$(getCatCmd "${input_fnames[$i]}")
done
for i in $(seq 0 $((n_update-1)));do
    cats["${update_fnames[$i]}"]=$(getCatCmd "${update_fnames[$i]}")
done

# get release, either from command line or from input file, in case of update
if [[ $n_update -eq 0 ]];then
    # if no release specified, set it to "1"
    if [[ -z "$release" ]];then
	release="1"
    fi
else
    # read previous release from input[0]
    release_col=$(getColNum "${input_fnames[0]}" "RELEASE" ${cats[${input_fnames[0]}]})
    if [[ -z "$release_col" ]];then
	echo "ERROR: no RELEASE column found in the input file ${input_fnames[0]}" 1>&2
	exit 1
    fi
    created_col=$(getColNum "${input_fnames[0]}" "CREATED" ${cats[${input_fnames[0]}]})
    if [[ -z "$created_col" ]];then
	echo "ERROR: no CREATED column found in the input file ${input_fnames[0]}" 1>&2
	exit 1
    fi
    if [[ -z "$release" ]];then
	release=$(${cats[${input_fnames[0]}]} "${input_fnames[0]}"|cut -f ${release_col}|head -n 3|tail -n 1)
	release=$((release+1))
    fi
fi

outfile="${out_dir}/phenotypes_r${release}.txt.gz"
missfile="${out_dir}/stats_r${release}.txt.gz"

if [[ -f "$outfile" ]];then
    echo "ERROR: output file $outfile already exists" 1>&2
    exit 1
fi
logfile="${out_dir}/phenotypes_r${release}.log"

: > "$logfile"
date "+%F %H-%M-%S"|tee -a "$logfile"
echo "Current dir: ${PWD}"|tee -a "$logfile"
echo "Command line: $scriptname ${args[@]}"|tee -a "$logfile"
echo ""|tee -a "$logfile"
echo "Output release: $release"|tee -a "$logfile"
echo "Output file: $outfile"|tee -a "$logfile"
echo "INFO: join mode: $mode"|tee -a "$logfile"
echo ""|tee -a "$logfile"
#----------------------------------------------------

# checking if all rows have the same number of fields
for i in $(seq 0 $((n_input-1)));do
    checkFields "${input_fnames[$i]}" "${cats[${input_fnames[$i]}]}" "$logfile"
    checkDuplicatesHeader "${input_fnames[$i]}" "${cats[${input_fnames[$i]}]}" "$logfile"
    checkRow "${input_fnames[$i]}" 1 "${cats[${input_fnames[$i]}]}" "$logfile"
    echo ""|tee -a "$logfile"
done
for i in $(seq 0 $((n_update-1)));do
    checkFields "${update_fnames[$i]}" "${cats[${update_fnames[$i]}]}" "$logfile"
    checkDuplicatesHeader "${update_fnames[$i]}" "${cats[${update_fnames[$i]}]}" "$logfile"
    checkRow "${update_fnames[$i]}" 1 "${cats[${update_fnames[$i]}]}" "$logfile"
    echo ""|tee -a "$logfile"
done

# check classes row if we're updating
if [[ $n_update -gt 0 ]];then
    checkRow "${input_fnames[0]}" 2 "${cats[${input_fnames[0]}]}" "$logfile"
fi
#----------------------------------------------------

# get ID column 
for i in $(seq 0 $((n_input-1)));do
    x=$(getColNum "${input_fnames[$i]}" "$id_field" ${cats[${input_fnames[$i]}]})
    if [[ -z $x ]];then
	echo "ERROR: no \"$id_field\" found in ${input_fnames[$i]}"|tee -a "$logfile"
	exit 1
    fi    
    input_ID_column+=($x)
    input_nrows+=($(${cats[${input_fnames[$i]}]} "${input_fnames[$i]}"|wc -l))
    input_ncols+=($(${cats[${input_fnames[$i]}]} "${input_fnames[$i]}"|head -n 1|tr '\t' '\n'|wc -l))
    checkColumn "${input_fnames[$i]}" $input_ID_column ${cats[${input_fnames[$i]}]} "$logfile"
    checkDuplicatesColumn "${input_fnames[$i]}" $input_ID_column ${cats[${input_fnames[$i]}]} "$logfile"
    echo ""  |tee -a "$logfile"
done
for i in $(seq 0 $((n_update-1)));do
    x=$(getColNum "${update_fnames[$i]}" "$id_field" ${cats[${update_fnames[$i]}]})
    if [[ -z $x ]];then
	echo "ERROR: no \"$id_field\" found in ${update_fnames[$i]}"|tee -a "$logfile"
	exit 1
    fi    
    update_ID_column+=($x)
    update_nrows+=($(${cats[${update_fnames[$i]}]} "${update_fnames[$i]}"|wc -l))
    update_ncols+=($(${cats[${update_fnames[$i]}]} "${update_fnames[$i]}"|head -n 1|tr '\t' '\n'|wc -l))
    checkColumn "${update_fnames[$i]}" $update_ID_column ${cats[${update_fnames[$i]}]} "$logfile"
    checkDuplicatesColumn "${update_fnames[$i]}" $update_ID_column ${cats[${update_fnames[$i]}]} "$logfile"
    echo ""  |tee -a "$logfile"
done
#----------------------------------------------------

# output info
for i in $(seq 0 $((n_input-1)));do
    echo "INFO: input file: ${input_fnames[$i]}"|tee -a "$logfile"
    echo "INFO: ID column index: ${input_ID_column[$i]}"|tee -a "$logfile"
    echo "INFO: rows: ${input_nrows[$i]}"|tee -a "$logfile"
    echo "INFO: columns: ${input_ncols[$i]}"|tee -a "$logfile"
    echo ""|tee -a "$logfile"
done

for i in $(seq 0 $((n_update-1)));do
    echo "INFO: update file: ${update_fnames[$i]}"|tee -a "$logfile"
    echo "INFO: ID column index: ${update_ID_column[$i]}"|tee -a "$logfile"
    echo "INFO: rows: ${update_nrows[$i]}"|tee -a "$logfile"
    echo "INFO: columns: ${update_ncols[$i]}"|tee -a "$logfile"
    echo ""|tee -a "$logfile"
done
#----------------------------------------------------

#
# check if column names (excepth for ID field) in input files are disjoint
#
echo -n "Checking if column names in input files are disjoint ... "|tee -a "$logfile"
if [[ $n_input -gt 1 ]];then
    for i in $(seq 0 $((n_input-1)));do
	for j in $(seq $((i+1)) $((n_input-1)));do
	    x=$(cat <("${cats[${input_fnames[$i]}]}" "${input_fnames[$i]}"|head -n 1|cut --complement -f ${input_ID_column[$i]}) <("${cats[${input_fnames[$j]}]}" "${input_fnames[$j]}"|head -n 1|cut --complement -f ${input_ID_column[$j]})|sort|uniq -d|wc -l)
	    if [[ $x -ne 0 ]];then
		echo "ERROR: input files ${input_fnames[$i]} and ${input_fnames[$j]} have columns in common"|tee -a "$logfile"
		exit 1
	    fi	    
	done
    done
fi
echo "OK"|tee -a "$logfile"
#----------------------------------------------------

#
# check if column names (excepth for ID field) in update files are disjoint
#
echo -n "Checking if column names in update files are disjoint ... "|tee -a "$logfile"
if [[ $n_update -gt 1 ]];then
    for i in $(seq 0 $((n_update-1)));do
	for j in $(seq $((i+1)) $((n_update-1)));do
	    x=$(cat <("${cats[${update_fnames[$i]}]}" "${update_fnames[$i]}"|head -n 1|cut --complement -f ${update_ID_column[$i]}) <("${cats[${update_fnames[$j]}]}" "${update_fnames[$j]}"|head -n 1|cut --complement -f ${update_ID_column[$j]})|sort|uniq -d|wc -l)
	    if [[ $x -ne 0 ]];then
		echo "ERROR: update files ${update_fnames[$i]} and ${update_fnames[$j]} have columns in common"|tee -a "$logfile"
		exit 1
	    fi	    
	done
    done
fi
echo "OK"|tee -a "$logfile"
echo ""|tee -a "$logfile"

#-------------------------------------- OUTPUT -------------------------------------------------

if [[ $n_update -eq 0 ]];then
    # JUST MERGING INPUT FILES IF MORE THAN ONE SPECIFIED
    # OR
    # EXCLUDE IDS IN -x LIST IF ONLY ONE INPUT FILE SPECIFIED

    if [[ $n_input -eq 1 ]];then
	paste <(${cats["${input_fnames[0]}"]} "${input_fnames[0]}"|head -n 1) <(echo RELEASE CREATED|tr ' ' '\t')|gzip - -c > "${outfile}"
    	for c in $("${cats[${input_fnames[0]}]}" "${input_fnames[0]}"|head -n 1|cut --complement -f ${input_ID_column[0]});do
    	    column_class[$c]="0"
    	done
	column_class["$id_field"]="NA"
	for c in $(${cats["${input_fnames[0]}"]} "${input_fnames[0]}"|head -n 1);do
	    class_array+=(${column_class[$c]})
	done
	class_array+=("NA" "NA")
	cline=$(join_by , "${class_array[@]}")
	eval "cat <(echo $cline|tr ',' '\t')|gzip - -c >> ${outfile}"

	${cats["${input_fnames[0]}"]} "${input_fnames[0]}" | perl -snle 'BEGIN{%h=();if (length($f)!=0){open(fh,"<",$f);while(<fh>){chomp;$h{$_}=1;}close(fh);}}{@a=split(/\t/);print $_ if !defined($h{$a[$c-1]});}' -- -f="$exclude_list" -c="${input_ID_column[0]}" | gzip - -c >> "${outfile}"
    else	
	# column classes of input columns are based on source input files
	for i in $(seq 0 $((n_input-1)));do
    	    for c in $("${cats[${input_fnames[$i]}]}" "${input_fnames[$i]}"|head -n 1|cut --complement -f ${input_ID_column[$i]});do
    		column_class[$c]=$i
    	    done
	done
	column_class["$id_field"]="NA"

	# OUTER JOIN on IDs
	tmpfile1=$(mktemp -p "$out_dir" temp_1_join_XXXXXXXX)
	if [[ $? -ne 0 ]];then
	    echo "ERROR: could not create tmpfile1 "|tee -a "$logfile"
	    exit 1
	fi

	fmt="1.${input_ID_column[0]},2.${input_ID_column[1]}"
	for i in $(seq 2 ${input_ncols[0]});do
	    fmt=$fmt",1.$i"
	done
	for i in $(seq 2 ${input_ncols[1]});do
	    fmt=$fmt",2.$i"
	done
	
	join_cmd="join --header -t$'\t' -1 ${input_ID_column[0]} -2 ${input_ID_column[1]} -a 1 -a 2 -e NA -o $fmt <(cat <(${cats[${input_fnames[0]}]} ${input_fnames[0]}|head -n 1) <(${cats[${input_fnames[0]}]} ${input_fnames[0]}|tail -n +2|sort -k${input_ID_column[0]},${input_ID_column[0]})) <(cat <(${cats[${input_fnames[1]}]} ${input_fnames[1]}|head -n 1) <(${cats[${input_fnames[1]}]} ${input_fnames[1]}|tail -n +2|sort -k${input_ID_column[1]},${input_ID_column[1]})) | gawk -v FS='\t' -v OFS='\t' '{if (\$2==\"NA\"){\$2=\$1;}print \$0;}'| cut -f 2-"

	n=${input_ncols[0]}
	m=${input_ncols[1]}
	n=$((n+m-1))
	for i in $(seq 2 $((n_input-1)));do
	    fmt="1.1,2.${input_ID_column[$i]}"
	    for j in $(seq 2 $n);do
		fmt=$fmt",1.$j"
	    done	
	    for j in $(seq 2 ${input_ncols[$i]});do
		fmt=$fmt",2.$j"
	    done
	    
	    join_cmd=$join_cmd"|join --header -t$'\t' -1 1 -2 ${input_ID_column[$i]} -a 1 -a 2 -e NA -o $fmt - <(cat <(${cats[${input_fnames[$i]}]} ${input_fnames[$i]}|head -n 1) <(${cats[${input_fnames[$i]}]} ${input_fnames[$i]}|tail -n +2|sort -k${input_ID_column[$i]},${input_ID_column[$i]})) | gawk -v FS='\t' -v OFS='\t' '{if (\$2==\"NA\"){\$2=\$1;}print \$0;}'| cut -f 2-"
	    m=${input_ncols[$i]}
	    n=$((n+m-1))
	done
	eval "$join_cmd > $tmpfile1"

	# adding class line
	paste <(head -n 1 "$tmpfile1") <(echo RELEASE CREATED|tr ' ' '\t')|gzip - -c > "${outfile}"
	for c in $(head -n 1 "$tmpfile1");do
	    class_array+=(${column_class[$c]})
	done
	class_array+=("NA" "NA")
	cline=$(join_by , "${class_array[@]}")
	eval "cat <(echo $cline|tr ',' '\t')|gzip - -c >> ${outfile}"
	# adding body, excluding IDs from exclude list
	x=$(cat $tmpfile1|wc -l)
	x=$((x-1))
	paste <(tail -n +2 "$tmpfile1") <(yes $release $datestr | tr ' ' '\t' | head -n $x) | perl -snle 'BEGIN{%h=();if (length($f)!=0){open(fh,"<",$f);while(<fh>){chomp;$h{$_}=1;}close(fh);}}{@a=split(/\t/);print $_ if !defined($h{$a[$c-1]});}' -- -f="$exclude_list" -c=1|gzip - -c >> "${outfile}"
	if [[ $debug == "NO" ]];then
	    rm -f "$tmpfile1"
	fi
    fi
    
    echo "Done"|tee -a "$logfile"
    date "+%F %H-%M-%S" |tee -a "$logfile"
else
    # updating the first input file using update files
    
    echo "Updating input ... "|tee -a "$logfile"

    #----------------------------------------------------------------
    # merging all update files using OUTER JOIN and saving the result in a temporary file, with classes
    tmpfile0=$(mktemp -p "$out_dir" temp_0_join_XXXXXXXX)
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not create tmpfile0 "|tee -a "$logfile"
	exit 1
    fi
    tmpfile1=$(mktemp -p "$out_dir" temp_1_update_XXXXXXXX)
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not create tmpfile1 "|tee -a "$logfile"
	exit 1
    fi

    echo -n "Merging update files ... "|tee -a "$logfile"

    # column classes of update files
    for i in $(seq 0 $((n_update-1)));do
    	for c in $("${cats[${update_fnames[$i]}]}" "${update_fnames[$i]}"|head -n 1|cut --complement -f ${update_ID_column[$i]});do
    	    column_class[$c]=$i
    	done
    done
    column_class["$id_field"]="NA"

    if [[ $n_update -gt 1 ]];then
	fmt="1.${update_ID_column[0]},2.${update_ID_column[1]}"
	for i in $(seq 2 ${update_ncols[0]});do
	    fmt=$fmt",1.$i"
	done
	for i in $(seq 2 ${update_ncols[1]});do
	    fmt=$fmt",2.$i"
	done
	
	join_cmd="join --header -t$'\t' -1 ${update_ID_column[0]} -2 ${update_ID_column[1]} -a 1 -a 2 -e NA -o $fmt <(cat <(${cats[${update_fnames[0]}]} ${update_fnames[0]}|head -n 1) <(${cats[${update_fnames[0]}]} ${update_fnames[0]}|tail -n +2|sort -k${update_ID_column[0]},${update_ID_column[0]})) <(cat <(${cats[${update_fnames[1]}]} ${update_fnames[1]}|head -n 1) <(${cats[${update_fnames[1]}]} ${update_fnames[1]}|tail -n +2|sort -k${update_ID_column[1]},${update_ID_column[1]})) | gawk -v FS='\t' -v OFS='\t' '{if (\$2==\"NA\"){\$2=\$1;}print \$0;}'| cut -f 2-"
	n=${update_ncols[0]}
	m=${update_ncols[1]}
	n=$((n+m-1))
	for i in $(seq 2 $((n_update-1)));do
	    fmt="1.1,2.${update_ID_column[$i]}"
	    for j in $(seq 2 $n);do
		fmt=$fmt",1.$j"
	    done	
	    for j in $(seq 2 ${update_ncols[$i]});do
		fmt=$fmt",2.$j"
	    done
	    
	    join_cmd=$join_cmd"|join --header -t$'\t' -1 1 -2 ${update_ID_column[$i]} -a 1 -a 2 -e NA -o $fmt - <(cat <(${cats[${update_fnames[$i]}]} ${update_fnames[$i]}|head -n 1) <(${cats[${update_fnames[$i]}]} ${update_fnames[$i]}|tail -n +2|sort -k${update_ID_column[$i]},${update_ID_column[$i]})) | gawk -v FS='\t' -v OFS='\t' '{if (\$2==\"NA\"){\$2=\$1;}print \$0;}'| cut -f 2-"
	    m=${update_ncols[$i]}
	    n=$((n+m-1))
	done
	eval "$join_cmd > $tmpfile0"
    else
	cp "${update_fnames[0]}" "$tmpfile0"
    fi
    
    head -n 1 "$tmpfile0" > "$tmpfile1"
    for c in $(head -n 1 "$tmpfile0");do
    	class_array+=(${column_class[$c]})
    done
    cline=$(join_by , "${class_array[@]}")
    echo $cline|tr ',' '\t' >> "$tmpfile1"
    tail -n +2 "$tmpfile0" >> "$tmpfile1"    
    if [[ $debug == "NO" ]];then
	rm -f "$tmpfile0"
    fi
	
    echo "Done"|tee -a "$logfile"
    #----------------------------------------------------------------
    # input fields (except for ID column) and their classes, from the input file
    declare -a input_fields
    declare -a input_classes
    declare -A input_field2class
    for c in $("${cats[${input_fnames[0]}]}" "${input_fnames[0]}"|head -n 1|cut --complement -f ${input_ID_column[0]},${release_col},${created_col});do
	input_fields+=($c)
    done
    for c in $("${cats[${input_fnames[0]}]}" "${input_fnames[0]}"|head -n 2|tail -n 1|cut --complement -f ${input_ID_column[0]},${release_col},${created_col});do
	input_classes+=($c)
    done
    for (( i=0; i<"${#input_fields[@]}"; i++ )); do input_field2class[${input_fields[$i]}]=${input_classes[$i]};done
    new_update_ID=$(getColNum "$tmpfile1" "$id_field" "cat")
    if [[ -z $new_update_ID ]];then
	echo "ERROR: no \"$id_field\" found in $tmpfile1"|tee -a "$logfile"
	exit 1
    fi    
    # update fields (except for ID column) and their classes, from the merged update file
    declare -a update_fields
    declare -a update_classes
    declare -A update_field2class
    for c in $(head -n 1 "$tmpfile1"|cut --complement -f $new_update_ID);do update_fields+=($c);done
    for c in $(head -n 2 "$tmpfile1"|tail -n 1|cut --complement -f $new_update_ID);do update_classes+=($c);done
    for (( i=0; i<"${#update_fields[@]}"; i++ )); do update_field2class[${update_fields[$i]}]=${update_classes[$i]};done
    
    # input classes intersecting with update fields: all input fields from these will be excluded from input before merging
    declare -A intersecting_classes
    # new fields in update that will be added
    declare -A new_fields
    for c in $(head -n 1 "$tmpfile1"|cut --complement -f $new_update_ID);do
	status=${input_field2class[$c]}
	if [[ ! -z "$status" ]];then
	    x=${input_field2class[$c]}
	    intersecting_classes[$x]=1
	else
	    new_fields[$c]=1
	fi
    done
    
    # all column numbers from input intersecting_classes
    declare -a input_columns_to_exclude
    for (( i=0; i<"${#input_fields[@]}"; i++ )); do
	c=${input_fields[$i]}
	x=${input_field2class[$c]}
	status=${intersecting_classes[$x]}
	if [[ ! -z "$status" ]];then
	    n=$(getColNum ${input_fnames[0]} $c ${cats[${input_fnames[0]}]})
	    input_columns_to_exclude+=($n)
	fi
    done
    input_columns_to_exclude+=($created_col)
    input_columns_to_exclude+=($release_col)

    # input without excluded columns, without classes
    tmpfile2=$(mktemp -p "$out_dir" temp_2_input_XXXXXXXX)
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not create tmpfile2 "|tee -a "$logfile"
	if [[ $debug == "NO" ]];then
	    rm -f "$tmpfile1"
	fi
	exit 1
    fi
    
    echo -n "Removing columns from input ... "|tee -a "$logfile"
    "${cats[${input_fnames[0]}]}" "${input_fnames[0]}"|cut -f $(join_by , "${input_columns_to_exclude[@]}") --complement| gawk 'BEGIN{FS=OFS="\t";}NR!=2{print $0;}' > "$tmpfile2"
    echo "Done"|tee -a "$logfile"
    new_input_ID=$(getColNum "$tmpfile2" "$id_field" "cat")
    if [[ -z $new_input_ID ]];then
	echo "ERROR: no \"$id_field\" found in $tmpfile2"|tee -a "$logfile"
	exit 1
    fi    

    common_IDs=$(join -1 1 -2 1 <(cut -f $new_update_ID "$tmpfile1"|tail -n +3|sort) <(cut -f $new_input_ID "$tmpfile2"|tail -n +2|sort)|wc -l)
    input_only_IDs=$(join -a 2 -e NULL -o 1.1,2.1 -1 1 -2 1 <(cut -f $new_update_ID "$tmpfile1"| tail -n +3|sort) <(cut -f $new_input_ID $tmpfile2|tail -n +2|sort)|grep NULL|wc -l)
    update_only_IDs=$(join -a 1 -e NULL -o 1.1,2.1 -1 1 -2 1 <(cut -f $new_update_ID "$tmpfile1"| tail -n +3|sort) <(cut -f $new_input_ID $tmpfile2|tail -n +2|sort)|grep NULL|wc -l)
    echo ""|tee -a "$logfile"
    echo "INFO: common IDs between input and update: $common_IDs"|tee -a "$logfile"
    echo "INFO: IDs only in input: $input_only_IDs (not included in output)"|tee -a "$logfile"
    str=""
    if [[ $mode == "inner" ]];then
	str="(not included in output, mode=$mode)"
    else
	str="(included in output, mode=$mode)"
    fi
    echo "INFO: IDs only in update: ${update_only_IDs} $str"|tee -a "$logfile"
    echo ""|tee -a "$logfile"

    fmtstr="1.""${new_input_ID}"",2.""${new_update_ID}"
    input_ncol=$(head -n 1 "$tmpfile2"|tr '\t' '\n'|wc -l)
    update_ncol=$(head -n 1 "$tmpfile1"|tr '\t' '\n'|wc -l)
    for (( i=1; i<=$input_ncol; i++ )); do
	if [[ $i -ne $new_input_ID ]];then
	    fmtstr=$fmtstr",1.$i"
	fi
    done
    for (( i=1; i<=$update_ncol; i++ )); do
	if [[ $i -ne $new_update_ID ]];then
	    fmtstr=$fmtstr",2.$i"
	fi
    done
    
    # join input and merged update, without classes
    tmpfile3=$(mktemp -p "$out_dir" temp_3_joined_XXXXXXXX)
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not create tmpfile3 "|tee -a "$logfile"
	if [[ $debug == "NO" ]];then
	    rm -f "$tmpfile1"
	    rm -f "$tmpfile2"
	fi
	exit 1
    fi
    
    echo -n "Joining ... "|tee -a "$logfile"
    join --header -1 $new_input_ID -2 $new_update_ID -a 1 -a 2 -t $'\t' -e NA -o $fmtstr "$tmpfile2" <(gawk 'BEGIN{FS=OFS="\t";}NR!=2{print $0;}' "$tmpfile1")|gawk -v m=$mode 'BEGIN{FS=OFS="\t";}{if (NR==1){print $0;}else{if ($2!="NA"){if (m=="inner"){if ($1!="NA"){print $0;}} if(m=="right"){if($1=="NA"){$1=$2;}print $0;} } }}'| cut -f 2 --complement > "$tmpfile3"
    echo "Done"|tee -a "$logfile"
    joined_ID=$(getColNum $tmpfile3 "$id_field" "cat")
    if [[ -z $joined_ID ]];then
	echo "ERROR: no \"$id_field\" found in $tmpfile3"|tee -a "$logfile"
	exit 1
    fi    

    # max update class
    maxc=0
    for c in "${update_field2class[@]}";do
	if [[ $c -gt $maxc ]];then
	    maxc=$c
	fi
    done
    declare -A new_classes
    declare -A tmp_ar
    cur_add=1
    for c in $(head -n 1 "$tmpfile3"|cut -f $joined_ID --complement);do
	s=${update_field2class[$c]}
	if [[ ! -z "$s" ]];then
	    new_classes[$c]=${update_field2class[$c]}
	else
	    t=${input_field2class[$c]}
	    if [[ ! -z "$t" ]];then
		d=${input_field2class[$c]}
		u=${tmp_ar[$d]}
		if [[ ! -z "$u" ]];then
		    new_classes[$c]=${tmp_ar[$d]}
		else
		    tmp_ar[$d]=$((maxc+cur_add))
		    new_classes[$c]=${tmp_ar[$d]}
		    cur_add=$((cur_add+1))
		fi
	    else
		echo "ERROR: could not find class for $c"|tee -a "$logfile"
		exit 1
	    fi
	fi
    done

    # FINAL OUTPUT
    # excluding IDs in exclude list
    echo -n "Creating output file ... "|tee -a "$logfile"
    command1="paste <(head -n 1 $tmpfile3) <(echo RELEASE CREATED|tr ' ' '\t')"

    declare -a temp
    for c in $(head -n 1 "$tmpfile3");do
	if [[ $c == "$id_field" ]];then
	    temp+=("NA")
	else
	    temp+=(${new_classes[$c]})
	fi
    done
    temp+=("NA" "NA")
    header2=$(join_by , "${temp[@]}")

    command2="paste <(tail -n +2 $tmpfile3)"
    x=$(tail -n +2 "$tmpfile3"| wc -l)
    command2=${command2}" <(yes $release $datestr| tr ' ' '\t'| head -n $x)"

    eval "cat <($command1) <(echo $header2|tr ',' '\t') <($command2)" | perl -snle 'BEGIN{%h=();if (length($f)!=0){open(fh,"<",$f);while(<fh>){chomp;$h{$_}=1;}close(fh);}}{@a=split(/\t/);print $_ if !defined($h{$a[$c-1]});}' -- -f="$exclude_list" -c="$joined_ID" | gzip - -c > "${outfile}"
    echo "Done"|tee -a "$logfile"

    if [[ $debug == "NO" ]];then
	rm -f "$tmpfile1"
	rm -f "$tmpfile2"
	rm -f "$tmpfile3"
    fi	
fi

# missingness counts
if [[ "$skipmiss" == "NO" ]];then
    echo "Creating missingness counts"|tee -a "$logfile"
    
    tmpfile_miss=$(mktemp -p "$out_dir" temp_4_miss_XXXXXXXX)    
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not create tmpfile_miss "|tee -a "$logfile"
	exit 1
    fi
    zcat "$outfile"|head -n 1 > "$tmpfile_miss"
    zcat "$outfile"|tail -n +3|parallel --citation --pipe -j64 -L10000 "${collectstats2}"|datamash -s -g 1 sum 2 >> "$tmpfile_miss"
    cat "$tmpfile_miss"| "${collectstats}" | gzip - -c > "$missfile"
    if [[ $debug == "NO" ]];then
	rm -f "$tmpfile_miss"
    fi
fi

echo "Done"|tee -a "$logfile"
date "+%F %H-%M-%S" |tee -a "$logfile"
    
exit 0

