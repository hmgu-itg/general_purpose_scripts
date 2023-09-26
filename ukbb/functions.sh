
# get total number of CPUs
function totalCPUs {
    echo $(grep ^processor /proc/cpuinfo|wc -l)
}

# get fraction of free CPUs
# a CPU is free if its idle time percentage is > 95%
function getFreeCPUs {
    local frac=1.0
    if [[ $# -ne 0 ]];then frac=$1;fi
    if (( $(echo "$frac > 1.0"|bc -l) ));then frac=1.0;fi
    totalfree=$(mpstat -P ALL 1 1 | sed 's/  */ /g' | awk 'BEGIN{c=0;}/Average:/ && $2 ~ /[0-9]/ && $12>95.0 {c=c+1;}END{print c;}')
    echo "${totalfree}*${frac}"| bc -l | awk '{print int($1);}'
}

# create temp files
# upon success the input associative array contains names of the created files as values
function createTempFiles {
    local dirname=$1
    # associative array, values: templates for mktemp
    local -n local_fnames=$2

    declare -a temp
    g=0
    for k in "${!local_fnames[@]}";do
	t="${local_fnames[$k]}"
	f=$(mktemp -p "$dirname" "$t")
	if [[ $? -ne 0 ]];then
	    echo "ERROR: could not create temporary file" >&2
	    g=1
	    break
	else
	    #echo "INFO: created temporary file $f" >&2
	    local_fnames[$k]="$f"
	    temp+=("$f")
	fi
    done
    
    if [[ "$g" -eq 1 ]];then
	for f in "${temp[@]}";do
	    echo "INFO: removing $f" >&2
	    rm "$f"
	done
	return 1
    else
	return 0
    fi
}

# return an associative array with column names as keys and 1-based column numbers as values, input is tab separated
function getColNumbers {
    local fname=$1
    local cmd=$2
    local -n local_arname=$3

    local_arname=()

    while read i x;do
	local_arname[$x]=$i
    done < <(eval "$cmd $fname"|head -n 1|perl -lne '$,=" ";@a=split(/\t/,$_,-1);for ($i=0;$i<scalar(@a);$i++){print $i+1,$a[$i];}')
}

# return an associative array with column names as keys and 1-based column numbers as values, field separator is provided
function getColNumbersSep {
    local fname=$1
    local cmd=$2
    local sep=$3
    local -n local_arname=$4

    local_arname=()

    while read i x;do
	local_arname[$x]=$i
    done < <(eval "$cmd $fname"|head -n 1|perl -slane '$,=" ";@a=split(/$sep/,$_,-1);for ($i=0;$i<scalar(@a);$i++){print $i+1,$a[$i];}' -- -sep="$sep")
}

# tab separated input
# return an array with field names
function getColNames {
    local fname=$1
    local cmd=$2
    local -n arname=$3

    arname=()

    # without f.eid, CREATED, RELEASE
    while read x;do
	arname+=($x)
    done < <(eval "$cmd $fname"|head -n 1|perl -lne '$,=" ";@a=split(/\t/,$_,-1);%H=();for ($i=0;$i<scalar(@a);$i++){if ($a[$i]=~/^f\.(\d+)\.\d+\.\d+$/){$H{$1}=1;}}foreach $k (sort {$a <=> $b} keys %H){print $k;}')
}

# for a given field, return an associative array "column index" --> "column name" for columns with names matching the field
# 12345 matches f.12345.0.1
function getCols {
    local fname=$1
    local cmd=$2
    local field=$3
    local -n local_arname=$4

    local_arname=()
    
    while read i x;do
	local_arname[$i]=$x
    done < <(eval "$cmd $fname"|head -n 1|perl -slne '$,=" ";@a=split(/\t/,$_,-1);for ($i=0;$i<scalar(@a);$i++){if ($a[$i]=~/^f\.$f\.\d+\.\d+$/){print $i+1,$a[$i];}}' -- -f=$field)
}

# read variable from a tab separated file
# 1st field: key
# 2nd field: value
#
# arguments: filename,key,variable name
function readValue {
    local fname=$1
    local key=$2
    local -n local_name=$3

    local_name=$(cat "$fname"|grep -v "#"|grep -m 1 "^$key" "$fname"|cut -f 2)
}

# read associative array from a tab separated file
# 1st field: key
# 2nd field: encoded associative array k1:v1,k2:v2,k3:v3, ...
#
# arguments: filename,key,array name
function readAArray {
    local fname=$1
    local key=$2
    local -n local_arname=$3

    local_arname=()
    
    while read x z;do
	local_arname["$x"]="$z"
    done < <(cat "$fname"|grep -v "#"|grep -m 1 "^$key" "$fname"|cut -f 2|tr ',' '\n'|tr ':' ' ')
}

# print error message and exit if the first argument is empty
function exitIfEmpty {
    local x=$1
    local msg=$2
    if [[ -z "$x" ]];then
	echo $msg
	exit 1
    fi
}

# print error message and exit if the first argument is not a file
function exitIfNotFile {
    local fname=$1
    local msg=$2
    if [[ ! -f "$fname" ]];then
	echo $msg
	exit 1
    fi
}

# print error message and exit if the first argument (file or dir) does not exist
function exitUnlessExists {
    local fname=$1
    local msg=$2
    if [[ ! -e "$fname" ]];then
	echo $msg
	exit 1
    fi
}

# print error message and exit if the first argument is an existing dir or file
function exitIfExists {
    local fname=$1
    local msg=$2
    if [[ -f "$fname" || -d "$fname" ]];then
	echo $msg
	exit 1
    fi
}

# print error message and exit if the first argument is not a directory
function exitIfNotDir {
    local dname=$1
    local msg=$2
    if [[ ! -d "$dname" ]];then
	echo $msg
	exit 1
    fi
}

# check if necessary commands are present
function checkInstalledCommands {
    local cmds=("$@")
    s="0"
    for c in "${cmds[@]}";do
	command -v $c > /dev/null
	if [[ $? -ne 0 ]];then
	    echo "ERROR: command $c not found" 1>&2
	    s="1"
	else
	    echo "INFO: check if $c is installed: OK" 1>&2
	fi
    done
    echo "$s"    
}

# check if key is in array
function checkArray {
    local k=$1
    shift
    local ar=("$@")

    x="NO"
    for c in "${ar[@]}";do
	if [[ $c == $k ]];then
	    x="YES"
	    break
	fi
    done
    
    echo $x
}

# given a value and an array, get 1-based index (or 0 if not in array) of value in array
function getArrayIndex {
    local k=$1
    shift
    local ar=("$@")

    local x="0"
    for (( i=0; i<${#ar[@]}; i++ ));
    do
	if [[ "${ar[$i]}" -eq $k ]];then
	    x=$((i+1))
	    break
	fi
    done
    echo $x
}

# check for duplicate fields in a column
function checkDuplicatesColumn {
    local fname=$1
    local column=$2

    local cmd="cat"
    if [[ $# -ge 3 ]];then
	cmd=$3
    fi

    local logfile=""
    if [[ $# -ge 4 ]];then
	logfile=$4
    fi

    totcols=$($cmd "$fname"|awk -F'\t' '{print NF; exit}')
    if [[ $column -gt $totcols ]];then
	if [[ -z "$logfile" ]];then
	    echo "WARN: checkDuplicatesColumn: column specified ($column) is greater than the total number of columns ($totcols) in $fname" 1>&2
	else
	    echo "WARN: checkDuplicatesColumn: column specified ($column) is greater than the total number of columns ($totcols) in $fname"|tee -a "$logfile"
	fi
	return 1
    fi
    
    if [[ -z "$logfile" ]];then
	echo -n "Checking for duplicates in column $column in $fname ... " 1>&2
    else
	echo -n "Checking for duplicates in column $column in $fname ... "|tee -a "$logfile"
    fi

    local x=$($cmd "$fname"|cut -f $column|sort|uniq -d|wc -l)
    if [[ $x -eq 0 ]];then
	if [[ -z "$logfile" ]];then
	    echo "OK" 1>&2
	else
	    echo "OK" | tee -a "$logfile"
	fi    
    else
	if [[ -z "$logfile" ]];then
	    echo "ERROR: column $column in $fname contains duplicates" 1>&2
	else
	    echo "ERROR: column $column in $fname contains duplicates"|tee -a "$logfile"
	fi
	exit 1
    fi
}

# check for duplicate fields in the header row
function checkDuplicatesHeader {
    local fname=$1
    local cmd="cat"
    if [[ $# -ge 2 ]];then
	cmd=$2
    fi

    local logfile=""
    if [[ $# -ge 3 ]];then
	logfile=$3
    fi

    if [[ -z "$logfile" ]];then
	echo -n "Checking for duplicates in the header row in $fname ... " 1>&2
    else
	echo -n "Checking for duplicates in the header row in $fname ... "|tee -a "$logfile"
    fi

    local x=$($cmd "$fname"|head -n 1|tr '\t' '\n'|sort|uniq -d|wc -l)
    if [[ $x -eq 0 ]];then
	if [[ -z "$logfile" ]];then
	    echo "OK" 1>&2
	else
	    echo "OK" | tee -a "$logfile"
	fi    
    else
	if [[ -z "$logfile" ]];then
	    echo "ERROR: header in $fname contains duplicate field(s)" 1>&2
	else
	    echo "ERROR: header in $fname contains duplicate field(s)"|tee -a "$logfile"
	fi
	exit 1
    fi
}

# check for empty fields in a row
function checkRow {
    local fname=$1
    local row=$2

    local cmd="cat"
    if [[ $# -ge 3 ]];then
	cmd=$3
    fi

    local logfile=""
    if [[ $# -ge 4 ]];then
	logfile=$4
    fi

    totrows=$($cmd "$fname"|wc -l)
    if [[ $row -gt $totrows ]];then
	if [[ -z "$logfile" ]];then
	    echo "WARN: checkRow: row specified ($row) is greater than the total number of rows ($totrows) in $fname" 1>&2
	else
	    echo "WARN: checkRow: row specified ($row) is greater than the total number of rows ($totrows) in $fname"|tee -a "$logfile"
	fi
	return 1
    fi
        
    if [[ -z "$logfile" ]];then
	echo -n "Checking for empty fields in row $row in $fname ... " 1>&2
    else
	echo -n "Checking for empty fields in row $row in $fname ... "|tee -a "$logfile"
    fi

    local x=$($cmd "$fname"|head -n $row|tail -n 1|gawk 'BEGIN{FS="\t";c=0;}{for (i=1;i<=NF;i++){if ($i ~ /^ *$/){c=i;exit;}}}END{print c;}')
    if [[ $x -eq 0 ]];then
	if [[ -z "$logfile" ]];then
	    echo "OK" 1>&2
	else
	    echo "OK" | tee -a "$logfile"
	fi    
    else
	if [[ -z "$logfile" ]];then
	    echo "ERROR: row $row in $fname contains empty field(s)" 1>&2
	else
	    echo "ERROR: row $row in $fname contains empty field(s)"|tee -a "$logfile"
	fi
	exit 1
    fi
}

# check for empty fields in a column
function checkColumn {
    local fname=$1
    local column=$2

    local cmd="cat"
    if [[ $# -ge 3 ]];then
	cmd=$3
    fi

    local logfile=""
    if [[ $# -ge 4 ]];then
	logfile=$4
    fi

    totcols=$($cmd "$fname"|awk -F'\t' '{print NF; exit}')
    if [[ $column -gt $totcols ]];then
	if [[ -z "$logfile" ]];then
	    echo "WARN: checkColumn: column specified ($column) is greater than the total number of columns ($totcols) in $fname" 1>&2
	else
	    echo "WARN: checkColumn: column specified ($column) is greater than the total number of columns ($totcols) in $fname"|tee -a "$logfile"
	fi
	return 1
    fi
    
    if [[ -z "$logfile" ]];then
	echo -n "Checking for empty fields in column $column in $fname ... " 1>&2
    else
	echo -n "Checking for empty fields in column $column in $fname ... "|tee -a "$logfile"
    fi

    local x=$($cmd $fname|cut -f $column|tr '\n' '\t'|sed 's/\t$//'|gawk 'BEGIN{FS="\t";c=0;}{for (i=1;i<=NF;i++){if ($i ~ /^ *$/){c=i;exit;}}}END{print c;}')
    if [[ $x -eq 0 ]];then
	if [[ -z "$logfile" ]];then
	    echo "OK" 1>&2
	else
	    echo "OK" | tee -a "$logfile"
	fi    
    else
	if [[ -z "$logfile" ]];then
	    echo "ERROR: column $column in $fname contains empty field(s)" 1>&2
	else
	    echo "ERROR: column $column in $fname contains empty field(s)"|tee -a "$logfile"
	fi
	exit 1
    fi
}

# check if all rows in a TSV file have same # fields
function checkFields {
    local fname=$1

    local cmd="cat"
    if [[ $# -ge 2 ]];then
	cmd=$2
    fi

    local logfile=""
    if [[ $# -ge 3 ]];then
	logfile=$3
    fi

    if [[ -z "$logfile" ]];then
	echo -n "Checking # fields in $fname ... " 1>&2
    else
	echo -n "Checking # fields in $fname ... " | tee -a "$logfile"
    fi
    
    local x=$($cmd "$fname"|gawk 'BEGIN{FS="\t";}{print NF;}'|sort|uniq|wc -l)
    if [[ $x -eq 1 ]];then
	if [[ -z "$logfile" ]];then
	    echo "OK" 1>&2
	else
	    echo "OK" | tee -a "$logfile"
	fi    
    else
	if [[ -z "$logfile" ]];then
	    echo "ERROR: $fname contains rows with different number of fields" 1>&2
	else
	    echo "ERROR: $fname contains rows with different number of fields" | tee -a "$logfile"
	fi
	exit 1
    fi
}

# join an array
function join_by { local IFS="$1"; shift; echo "$*"; }

# get an array of column numbers for a specific column name (can be several columns with the same name)
function getColNums {
    local fname=$1
    local colname=$2
    local catcmd=$3
    local -n arname=$4

    arname=()
    while read i;do
	arname+=($i)
    done < <(eval "$catcmd $fname" | head -n 1 | perl -slne '@a=split(/\t/,$_,-1);for ($i=0;$i<scalar(@a);$i++){print $i+1 if ($n eq $a[$i]);}' -- -n="$colname")
}

# get column number for a specific column
function getColNum () {
    local fname=$1
    local colname=$2
    local cmd="cat"
    if [[ $# -ge 3 ]];then
	cmd=$3
    fi
    echo $(fgrep -w "$colname"  <($cmd "$fname"|head -n 1|tr '\t' '\n'|cat -n|sed 's/^  *//')|cut -f 1)
}

# get column number for a specific column in a file inside a tar.gz
function getTGZColNum () {
    local tgz_fname=$1
    local fname=$2
    local colname=$3
    echo $(fgrep -w "$colname"  <(tar -zxf "$tgz_fname" "$fname" -O|head -n 1|tr '\t' '\n'|cat -n|sed 's/^  *//')|cut -f 1)
}

# get cat/zcat, based on file name
function getCatCmd () {
    local fname=$1
    if [[ "$fname" =~ gz$ ]];then
	echo "zcat"
    else
	echo "cat"
    fi
}
