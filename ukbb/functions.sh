
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
    echo -n "Checking # fields in $fname ... " 1>&2
    local x=$(awk 'BEGIN{FS="\t";}{print NF;}' $fname| sort|uniq| wc -l)
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
    echo $(fgrep -w $colname  <(head -n 1 $fname | tr '\t' '\n'| cat -n | sed 's/^  *//') | cut -f 1)
}
