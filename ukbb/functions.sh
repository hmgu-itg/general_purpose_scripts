
# check if necessary commands are present
function checkInstalledCommands {
    local cmds=("gawk")
    s="0"
    for c in ${cmds[@]};do
	command -v $c > /dev/null
	if [[ $? -ne 0 ]];then
	    echo "ERROR: command $c not found" 1>&2
	    s="1"
	else
	    echo "INFO: check $c: OK" 1>&2
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

function getCatCmd () {
    local fname=$1
    if [[ "$fname" =~ gz$ ]];then
	echo "zcat"
    else
	echo "cat"
    fi
}
