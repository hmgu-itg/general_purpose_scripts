
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
    
    echo -n "Checking # fields in $fname ... " 1>&2
    local x=$($cmd $fname|awk 'BEGIN{FS="\t";}{print NF;}'| sort|uniq| wc -l)
    if [[ $x -eq 1 ]];then
	echo "OK" 1>&2
    else
	echo "ERROR: $fname contains rows with different number of fields" 1>&2
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
	cmd=$2
    fi
    echo $(fgrep -w $colname  <($cmd $fname|head -n 1|tr '\t' '\n'|cat -n |sed 's/^  *//')|cut -f 1)
}

function getCatCmd () {
    local fname=$1
    if [[ "$fname"~/gz$/ ]];then
	echo "zcat"
    else
	echo "cat"
    fi
}
