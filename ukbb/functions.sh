
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
    
    local x=$($cmd $fname|gawk 'BEGIN{FS="\t";}{print NF;}'| sort|uniq| wc -l)
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
    echo $(fgrep -w $colname  <($cmd $fname|head -n 1|tr '\t' '\n'|cat -n |sed 's/^  *//')|cut -f 1)
}

function getCatCmd () {
    local fname=$1
    if [[ "$fname" =~ gz$ ]];then
	echo "zcat"
    else
	echo "cat"
    fi
}
