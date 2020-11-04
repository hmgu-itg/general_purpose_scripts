
function parseCommandLine () {
    aggs=("$@")
    for i in "${arr[@]}";
    do
	echo "ARG: " $i
    done
}
