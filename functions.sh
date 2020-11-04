
function getvalue () {
    args=("$@")
    key=${args[-1]}
    unset args[-1]

    for ((i=0; i <= ${#args[@]}; i++)); do
	if [[ "${args[i]}" == "-$key" ]];then
	    echo "${args[i+1]}"
	    break
	fi
    done
}


function parseCommandLine () {
    declare -a A
    args=("$@")
    s=${args[-1]}
    unset args[-1]
    for ((i=0; i<${#s}; i++)); do
	x="${s:$i:1}"
	z=$(getvalue ${args[@]} $x)
	A+=("$z")	
    done
    
    echo ${A[@]}
}
