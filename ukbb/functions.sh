
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
# input: dirname
# input: associative array name with filenames templates
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
	    # echo "INFO: created temporary file $f" >&2
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

function createTempFilesN {
    local dirname=$1
    local n=$2
    local -n fnames=$3
    local pattern="temp_XXXXXXXX"
    if [[ $# -ge 4 ]];then
	pattern=$4
    fi
    local i

    fnames=()
    for (( i=0; i<$n; i++ ));do
	fnames[$i]=$pattern
    done

    createTempFiles $dirname $3
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

# print error message and exit if the first argument is not a directory
function exitIfDir {
    local name=$1
    local msg=$2
    if [[ -d "$name" ]];then
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
    for (( i=0; i<${#ar[@]}; i++ ));do
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

    local i=$(eval "$cmd $fname" | head -n 1 | perl -slne 'BEGIN{$k="";}{@a=split(/\t/,$_,-1);for ($i=0;$i<scalar(@a);$i++){if ($a[$i] eq $n){$k=$i+1;exit 0;}}}END{print $k;}' -- -n="$colname")
    echo $i
    
    # echo $(fgrep -w "$colname"  <($cmd "$fname"|head -n 1|tr '\t' '\n'|cat -n|sed 's/^  *//')|cut -f 1)
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
    if [[ "$fname" =~ \.gz$ ]];then
	echo "zcat"
    else
	echo "cat"
    fi
}

# get release from a merged file
function getRelease {
    local fname=$1
    local cmd=$2

    local rc=$(getColNum "$fname" "RELEASE" "$cmd")
    if [[ -z "$rc" ]];then
	echo ""
    else
	echo $("$cmd" "$fname"|head -n 3|tail -n 1|cut -f $rc)
    fi
}

# get release from a merged file
function getCreated {
    local fname=$1
    local cmd=$2

    local rc=$(getColNum "$fname" "CREATED" "$cmd")
    if [[ -z "$rc" ]];then
	echo ""
    else
	echo $("$cmd" "$fname"|head -n 3|tail -n 1|cut -f $rc)
    fi
}

# report samples missing in at least one of the files
function report_missing {
    # name of an array containing file names
    local -n fnames=$1
    local -n ct=$2
    local -n idcols=$3
    local logfile=$4
    local n="${#fnames[@]}"
    declare -a ar
    local i

    ar+=("MISSING" "ID")
    for i in "${!fnames[@]}";do
	echo "MISSING $((i+1)) ${fnames[$i]}" | tr ' ' '\t' >> "$logfile"
	ar+=("$((i+1))")
    done
    echo $(join_by " " "${ar[@]}") | tr ' ' '\t' >> "$logfile"

    ar=()
    for i in "${!fnames[@]}";do
	ar+=("<(${ct[${fnames[$i]}]} ${fnames[$i]} | tail -n +2 | cut -f ${idcols[$i]})")
    done
    local cmd="cat "$(join_by " " "${ar[@]}")" | sort | uniq"
    
    local fmt="1.1,2.1"
    local join_cmd="join -1 1 -2 1 -a 1 -a 2 -e NA -o $fmt <($cmd) <(${ct[${fnames[0]}]} ${fnames[0]} | tail -n +2 | cut -f ${idcols[0]} | sort -k1,1)"
    for i in $(seq 1 "$((n-1))");do
	fmt="1.1"
	for j in $(seq 2 $((i+1)));do
	    fmt=$fmt",1.$j"
	done
	fmt=$fmt",2.1"
	join_cmd="$join_cmd"" | join -1 1 -2 1 -a 1 -a 2 -e NA -o $fmt - <(${ct[${fnames[$i]}]} ${fnames[$i]} | tail -n +2 | cut -f ${idcols[$i]} | sort -k1,1)"
    done
    # echo "DEBUG: $join_cmd" | tee -a "$logfile"
    
    while read i;do
	echo "MISSING $i" | tr ' ' '\t' >> "$logfile"
    done < <(eval "$join_cmd" | grep "NA" | awk '{for (i=2;i<=NF;i++){if ($i=="NA"){$i="N";}else{$i="Y";}}print $0;}')
}

# no common fields in input files allowed
function merge_two_files {
    local fname1=$1
    local fname2=$2
    local tmpdir=$3
    local logfile=$4
    local -n ret=$5
    local i
    local x

    local cat1=$(getCatCmd "$fname1")
    local cat2=$(getCatCmd "$fname2")
    local idCol1=$(getColNum "$fname1" "f.eid" "$cat1")
    local idCol2=$(getColNum "$fname2" "f.eid" "$cat2")
    local ncols1=$("$cat1" "$fname1" | head -n 1 | tr '\t' '\n' | wc -l)
    local ncols2=$("$cat2" "$fname2" | head -n 1 | tr '\t' '\n' | wc -l)

    local tmpfile=$(mktemp -p "$tmpdir" merge_pair_XXXXXXXX)
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not create tmpfile" | tee -a "$logfile"
	ret=""
	return
    fi

    echo "INFO: merging"  | tee -a "$logfile"
    echo "INFO: file1: $fname1" | tee -a "$logfile"
    echo "INFO: file2: $fname2" | tee -a "$logfile"
    echo "INFO: output: $tmpfile" | tee -a "$logfile"
    echo ""  | tee -a "$logfile"
    
    local fmt="1.${idCol1},2.${idCol2}"
    for i in $(seq 1 ${ncols1});do
	if [[ "$i" -ne "${idCol1}" ]];then
	   fmt=$fmt",1.$i"
	fi
    done
    for i in $(seq 1 ${ncols2});do
	if [[ "$i" -ne "${idCol2}" ]];then
	    fmt=$fmt",2.$i"
	fi
    done

    echo -n "INFO: joining input files ... "  | tee -a "$logfile"
    
    local join_cmd="join --header -t$'\t' -1 ${idCol1} -2 ${idCol2} -a 1 -a 2 -e NA -o $fmt <(cat <(${cat1} ${fname1} | head -n 1) <(${cat1} ${fname1} | tail -n +2 | sort -T ${tmpdir} -t$'\t' -k${idCol1},${idCol1})) <(cat <(${cat2} ${fname2} | head -n 1) <(${cat2} ${fname2} | tail -n +2 | sort -T ${tmpdir} -t$'\t' -k${idCol2},${idCol2}))"
    eval "$join_cmd > $tmpfile"
    echo "done"  | tee -a "$logfile"

    x=$(tail -n +2  "$tmpfile" | awk 'BEGIN{FS="\t";c=0;}{if ($1!="NA" && $2=="NA"){c=c+1;}}END{print c;}')
    echo "INFO: samples in file1 only: $x"  | tee -a "$logfile"
    x=$(tail -n +2  "$tmpfile" | awk 'BEGIN{FS="\t";c=0;}{if ($1=="NA" && $2!="NA"){c=c+1;}}END{print c;}')
    echo "INFO: samples in file2 only: $x"  | tee -a "$logfile"
    x=$(tail -n +2  "$tmpfile" | awk 'BEGIN{FS="\t";c=0;}{if ($1!="NA" && $2!="NA"){c=c+1;}}END{print c;}')
    echo "INFO: samples in both file1 and file2: $x"  | tee -a "$logfile"
    x=$(head -n 1 "$tmpfile" | tr '\t' '\n' | wc -l )
    echo "INFO: total columns in joined file: $x"  | tee -a "$logfile"
    awk 'BEGIN{FS=OFS="\t";}{if ($2=="NA"){$2=$1;}print $0;}' "$tmpfile" | cut -f 2- | TMPDIR="${tmpdir}" sponge "$tmpfile"
    echo -e "\n---------------------------------------------------------\n" | tee -a "$logfile"
    
    ret="$tmpfile"
}

# split n into m parts
function split_int {
    local n=$1
    local m=$2
    local -n local_arname=$3
    local i
    local x
    local z

    local_arname=()

    if [[ $m -le 0 ]];then
	echo "ERROR: split_int: first argument must be positive"
	return -1
    fi
    
    if [[ $m -le 0 ]];then
	echo "ERROR: split_int: second argument must be positive"
	return -1
    fi
    
    x=$(( n / m ))
    z=$(( n - x * m))

    for (( i=0; i<${m}; i++ ));do
	if [[ $i -lt $z ]];then
	    local_arname+=($(( x + 1 )))
	else
	    local_arname+=($x)
	fi
    done
}

# input files are tab separated, with header
# performs outer join on the first column of both files
# the columns of the bigger of the two files is split in <parts> parts (default: no splitting)
function join_two_files {
    local logfile=$1
    local tmpdir=$2
    local infile1=$3
    local infile2=$4
    local outfile=$5
    local parts=1
    if [[ $# -ge 6 ]];then
	parts=$6
    fi

    echo "DEBUG: join_two_files: infile1: $infile1"  | tee -a "$logfile"
    echo "DEBUG: join_two_files: infile2: $infile2"  | tee -a "$logfile"
    echo "DEBUG: join_two_files: outfile: $outfile"  | tee -a "$logfile"
    echo "DEBUG: join_two_files: tmpdir: $tmpdir"  | tee -a "$logfile"
    echo "DEBUG: join_two_files: logfile: $logfile"  | tee -a "$logfile"
    echo "DEBUG: join_two_files: parts: $parts"  | tee -a "$logfile"
    
    local cat1=$(getCatCmd "$infile1")
    local cat2=$(getCatCmd "$infile2")
    local cat3=$(getCatCmd "$outfile")

    echo "DEBUG: join_two_files: $cat1 $infile1"  | tee -a "$logfile"
    echo "DEBUG: join_two_files: $cat2 $infile2"  | tee -a "$logfile"
    echo "DEBUG: join_two_files: $cat3 $outfile"  | tee -a "$logfile"
    
    local i j k tdir tf c1 c2
    declare -a temp
    declare -A tempf
    
    echo "DEBUG: join_two_files: FILE 1: $infile1" | tee -a "$logfile"
    echo "DEBUG: join_two_files: FILE 2: $infile2" | tee -a "$logfile"
    echo "DEBUG: join_two_files: OUTPUT: $outfile" | tee -a "$logfile"
    echo "DEBUG: join_two_files: parts: $parts" | tee -a "$logfile"

    if [[ ! -e "$infile1" ]];then
	if [[ -e "$infile2" ]];then
	    echo "DEBUG: join_two_files: FILE 1 doesn't exist, copying FILE 2 to $outfile" | tee -a "$logfile"
	    cp "$infile2" "$outfile"
	    return
	else
	    echo "DEBUG: join_two_files: input files don't exist; exit" | tee -a "$logfile"
	    return
	fi
    else
	if [[ ! -e "$infile2" ]];then
	    echo "DEBUG: join_two_files: FILE 2 doesn't exist, copying FILE 1 to $outfile" | tee -a "$logfile"
	    cp "$infile1" "$outfile"
	    return
	fi
    fi
    
    local nc1=$("$cat1" "$infile1" | head -n 1 | tr '\t' '\n' | wc -l)
    local nc2=$("$cat2" "$infile2" | head -n 1 | tr '\t' '\n' | wc -l)

    if [[ $parts -eq 1 ]];then
	local fmt="1.1,2.1"
	for i in $(seq 2 ${nc1});do
	    fmt=$fmt",1.$i"
	done
	for i in $(seq 2 ${nc2});do
	    fmt=$fmt",2.$i"
	done
	# echo "DEBUG: fmt: $fmt"
	# "$cat1" "$infile1"
	# "$cat2" "$infile2"
	# echo "DEBUG: joining, $infile1, $infile2, $outfile, $fmt"
	# keeping both ID columns in output
	if [[ "$cat3" == "cat" ]];then
	    join --header -t$'\t' -1 1 -2 1 -a 1 -a 2 -e "NA" -o "$fmt" <(cat <("$cat1" "$infile1" | head -n 1) <("$cat1" "$infile1" | tail -n +2 | sort -t$'\t' -k1,1)) <(cat <("$cat2" "$infile2" | head -n 1) <("$cat2" "$infile2" | tail -n +2 | sort -T ${tmpdir} -t$'\t' -k1,1)) > "$outfile"
	else
	    join --header -t$'\t' -1 1 -2 1 -a 1 -a 2 -e "NA" -o "$fmt" <(cat <("$cat1" "$infile1" | head -n 1) <("$cat1" "$infile1" | tail -n +2 | sort -t$'\t' -k1,1)) <(cat <("$cat2" "$infile2" | head -n 1) <("$cat2" "$infile2" | tail -n +2 | sort -T ${tmpdir} -t$'\t' -k1,1)) | gzip - > "$outfile"
	fi
    else
	tdir=$(dirname $outfile)
	createTempFilesN "$tdir" $(( parts + 1 )) tempf "temp_XXXXXXXX.gz"
	tf="${tempf[0]}"

	if [[ $nc1 -gt $nc2 ]];then
	    split_int $(( nc1 -1 )) $parts temp
	    echo "DEBUG: join_two_files: FILE 1: split "$(( nc1 -1 ))" columns into $parts chunks" | tee -a "$logfile"
	    echo "DEBUG: join_two_files: columns split into chunks:" | tee -a "$logfile"
	    for i in "${temp[@]}";do
		echo "DEBUG: $i" | tee -a "$logfile"
	    done
	    c1=2 # start column
	    k=1 # current temp output filename
	    for i in ${temp[@]};do
		c2=$(( c1 + i -1 ))
		$cat1 "$infile1" | cut -f 1,"${c1}"-"${c2}" | gzip - > "$tf"
		# join_two_files
		join_two_files "$logfile" "$tmpdir" "$tf" "$infile2" "${tempf[$k]}"
		# echo "DEBUG: i=$i"
		# zcat "${tempf[$k]}"
		c1=$(( c2 + 1 ))
		k=$(( k + 1 ))
	    done
	    # paste all temp files
	    # temp files have two ID columns
	    cp "${tempf[1]}" "$tf"
	    j=2
	    for (( i=2; i<=${parts}; i++ ));do
		k=$(( i - 2 ))
		k="${temp[$k]}"
		j=$(( j + k ))
		# echo "DEBUG: k=$k j=$j"
		paste <(zcat "$tf" | cut -f 1-"$j") <(zcat "${tempf[$i]}" | cut --complement -f 1,2) | gzip - > "${tempf[1]}"
		cp "${tempf[1]}" "$tf"
	    done	
	else
	    split_int $(( nc2 -1 )) $parts temp
	    echo "DEBUG: FILE 2: split "$(( nc2 -1 ))" columns into $parts chunks" | tee -a "$logfile"
	    echo "DEBUG: join_two_files: columns split into chunks:" | tee -a "$logfile"
	    for i in "${temp[@]}";do
		echo "DEBUG: $i" | tee -a "$logfile"
	    done
	    c1=2 # start column
	    k=1 # current temp output filename
	    for i in ${temp[@]};do
		c2=$(( c1 + i -1 ))
		$cat2 "$infile2" | cut -f 1,"${c1}"-"${c2}" | gzip - > "$tf"
		# join_two_files
		join_two_files "$logfile" "$tmpdir" "$infile1" "$tf" "${tempf[$k]}"
		# echo "DEBUG: i=$i"
		# zcat "${tempf[$k]}"
		c1=$(( c2 + 1 ))
		k=$(( k + 1 ))
	    done
	    # paste all temp files
	    # temp files have two ID columns
	    cp "${tempf[1]}" "$tf"
	    for (( i=2; i<=${parts}; i++ ));do
		k="${temp[$i]}"
		k=$(( k + 2 ))
		paste <(zcat "$tf") <(zcat "${tempf[$i]}" | cut --complement -f 1-"${k}") | gzip - > "${tempf[1]}"
		cp "${tempf[1]}" "$tf"
	    done	    
	fi
	
	if [[ "$cat3" == "cat" ]];then
	    gunzip -c "$tf" > "$outfile"
	else
	    cp "$tf" "$outfile"
	fi

	# delete temp files
	for k in "${tempf[@]}";do
	    if [[ -f "$k" ]];then
		rm "$k"
	    fi
	done	
    fi    
}

function update_file {
    local fname1=$1
    local fname2=$2
    local tmpdir=$3
    local logfile=$4
    local -n ret=$5
    local xn=$6
    local ko=$7
    local parts=1
    if [[ $# -ge 8 ]];then
	parts=$8
    fi
    
    local c i x

    declare -A colnames1
    declare -A colnames2
    declare -a common_cols
    declare -a tmp_ar
    declare -a exclude_cols
    declare -a common_colnum

    local cat1=$(getCatCmd "$fname1")
    local cat2=$(getCatCmd "$fname2")
    local idCol1=$(getColNum "$fname1" "f.eid" "$cat1")
    local idCol2=$(getColNum "$fname2" "f.eid" "$cat2")
    local ncols1=$("$cat1" "$fname1" | head -n 1 | tr '\t' '\n' | wc -l)
    local ncols2=$("$cat2" "$fname2" | head -n 1 | tr '\t' '\n' | wc -l)

    getColNumbers "$fname1" "$cat1" colnames1
    getColNumbers "$fname2" "$cat2" colnames2

    echo "INFO: updating"  | tee -a "$logfile"
    echo "INFO: file1: $fname1"  | tee -a "$logfile"
    echo "INFO: file2: $fname2"  | tee -a "$logfile"
    echo "INFO: logfile: $logfile"  | tee -a "$logfile"
    echo "INFO: tmpdir: $tmpdir"  | tee -a "$logfile"
    echo "INFO: exclude: $xn"  | tee -a "$logfile"
    echo "INFO: keep old values: $ko"  | tee -a "$logfile"
    echo "INFO: parts: $parts"  | tee -a "$logfile"
    echo ""  | tee -a "$logfile"
    
    for c in "${!colnames1[@]}";do
	if [[ -v "colnames2[$c]" && "$c" != "f.eid" ]];then
	    common_cols+=($c)
	fi
    done

    echo "INFO: common fields: ${#common_cols[@]}"  | tee -a "$logfile"
    echo ""  | tee -a "$logfile"
    
    local tmpfile=$(mktemp -p "$tmpdir" merge_XXXXXXXX)
    if [[ $? -ne 0 ]];then
	echo "ERROR: could not create tmpfile" | tee -a "$logfile"
	ret=""
	return
    fi

    echo "INFO: output: $tmpfile"  | tee -a "$logfile"
    echo ""  | tee -a "$logfile"

    if [[ $parts -eq 1 ]];then
	local fmt="1.${idCol1},2.${idCol2}"
	for i in $(seq 1 ${ncols1});do
	    if [[ "$i" -ne "${idCol1}" ]];then
		fmt=$fmt",1.$i"
	    fi
	done
	for i in $(seq 1 ${ncols2});do
	    if [[ "$i" -ne "${idCol2}" ]];then
		fmt=$fmt",2.$i"
	    fi
	done

	echo "INFO: joining input files"  | tee -a "$logfile"
	echo ""  | tee -a "$logfile"
	
	local join_cmd="join --header -t$'\t' -1 ${idCol1} -2 ${idCol2} -a 1 -a 2 -e NA -o $fmt <(cat <(${cat1} ${fname1} | head -n 1) <(${cat1} ${fname1} | tail -n +2 | sort -T ${tmpdir} -t$'\t' -k${idCol1},${idCol1})) <(cat <(${cat2} ${fname2} | head -n 1) <(${cat2} ${fname2} | tail -n +2 | sort -T ${tmpdir} -t$'\t' -k${idCol2},${idCol2}))"
	# echo "DEBUG: join command line: $join_cmd"
	eval "$join_cmd > $tmpfile"
	if [[ $? -eq 0 ]];then
	    echo "INFO: join OK"  | tee -a "$logfile"
	else
	    echo "ERROR: join failed"  | tee -a "$logfile"
	fi
    else
	join_two_files "$logfile" "$tmpdir" "$fname1" "$fname2" "$tmpfile" "$parts"	
    fi
    
    x=$(tail -n +2  "$tmpfile" | awk 'BEGIN{FS="\t";c=0;}{if ($1!="NA" && $2=="NA"){c=c+1;}}END{print c;}')
    echo "INFO: samples in file1 only: $x"  | tee -a "$logfile"
    x=$(tail -n +2  "$tmpfile" | awk 'BEGIN{FS="\t";c=0;}{if ($1=="NA" && $2!="NA"){c=c+1;}}END{print c;}')
    echo "INFO: samples in file2 only: $x"  | tee -a "$logfile"
    x=$(tail -n +2  "$tmpfile" | awk 'BEGIN{FS="\t";c=0;}{if ($1!="NA" && $2!="NA"){c=c+1;}}END{print c;}')
    echo "INFO: samples in both file1 and file2: $x"  | tee -a "$logfile"
    x=$(head -n 1 "$tmpfile" | tr '\t' '\n' | wc -l )
    echo "INFO: total columns in joined file: $x"  | tee -a "$logfile"
    
    if [[ $ko == "NO" ]];then
	awk 'BEGIN{FS=OFS="\t";}{if ($2=="NA"){$2=$1;}print $0;}' "$tmpfile" | cut -f 2- | TMPDIR="${tmpdir}" sponge "$tmpfile"
    else
	common_colnum=()
	for c in "${common_cols[@]}";do
	    getColNums "$tmpfile" "$c" "cat" tmp_ar
	    common_colnum+=("${tmp_ar[0]}")
	    common_colnum+=("${tmp_ar[1]}")
	done
	#echo "common_colnum:" $(join_by "," "${common_colnum[@]}")
	#cat "$tmpfile"
	cat "$tmpfile" | "${scriptdir}/keep_old_values.pl" $(join_by "," "${common_colnum[@]}") | awk 'BEGIN{FS=OFS="\t";}{if ($2=="NA"){$2=$1;}print $0;}' | cut -f 2- | TMPDIR="${tmpdir}" sponge "$tmpfile"
    fi
    echo -e "\n---------------------------------------------------------\n" | tee -a "$logfile"
    # tmpfile contains one ID column (1st), all columns from file1 (including CREATED, RELEASE), all columns from file2
    
    if [[ "${#common_cols[@]}" -eq 0 ]];then # no common colnames, remove RELEASE, CREATED columns
	cut --complement -f $(getColNum "$tmpfile" "RELEASE"),$(getColNum "$tmpfile" "CREATED") "$tmpfile" | TMPDIR="${tmpdir}" sponge "$tmpfile"
	if [[ ! -z "$xn" ]];then
	    echo "INFO: no common columns, nothing to save"  | tee -a "$logfile"
	fi
	ret="$tmpfile"
    else # there are common colnames
	common_colnum=(1)
	for c in "${common_cols[@]}";do
	    getColNums "$tmpfile" "$c" "cat" tmp_ar
	    common_colnum+=("${tmp_ar[0]}")
	    common_colnum+=("${tmp_ar[1]}")
	done	
	echo "INFO: samples with x -> NA" >> "$logfile"
	while read line;do # for fields with x -> NA, report "Y", otherwise "N"
	    echo "x2NA" "$line" | tr ' ' '\t' >> "$logfile"
	done < <(cut -f $(join_by "," "${common_colnum[@]}") "$tmpfile" | "${scriptdir}/sort_columns.pl" | perl -lne '@a=split(/\t/,$_,-1);$f=0;@b=($a[0]);if ($_=~/^f/){$f=1;for ($i=1;$i<scalar(@a);$i+=2){push @b,$a[$i];}}else{for ($i=1;$i<scalar(@a);$i+=2){if ($a[$i] ne "NA" && $a[$i+1] eq "NA"){$f=1;push @b,"Y";}else{push @b,"N";};}}print join(" ",@b) if ($f==1);')
	if [[ ! -z "$xn" ]];then
	    echo "INFO: saving common columns to $xn" | tee -a "$logfile"
	    cut -f $(join_by "," "${common_colnum[@]}") "$tmpfile" | "${scriptdir}/sort_columns.pl" | gzip - -c > "$xn"
	fi
	# report some stats comparing shared columns in both files
	"${scriptdir}/get_stats.pl" "$tmpfile" >> "$logfile"
	# remove RELEASE, CREATED columns
	exclude_cols=()
	exclude_cols+=($(getColNum "$tmpfile" "RELEASE"))
	exclude_cols+=($(getColNum "$tmpfile" "CREATED"))
	# remove common columns coming from the first input file
	for c in "${common_cols[@]}";do
	    getColNums "$tmpfile" "$c" "cat" tmp_ar
	    exclude_cols+=("${tmp_ar[0]}")
	done
	cut --complement -f $(join_by "," "${exclude_cols[@]}") "$tmpfile" | TMPDIR="${tmpdir}" sponge "$tmpfile"

	ret="$tmpfile"
    fi
}

# tab separated input
# check if a value occurs in a column
function find_value_in_column {
    local fname=$1
    local col=$2
    local val=$3

    while read x;do
	if [[ "$x" == "$val" ]];then
	    return 0
	fi
    done < <(cut -f "$col" "$fname")

    return 1
}
