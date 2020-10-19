#!/bin/bash

max_array_size=$(grep "^MaxArraySize=" /etc/slurm/slurm.conf | cut -d '=' -f 2)

function usage
{
    echo "Usage: $0"
    echo "       -i <file with commands> : required"
    echo "       -m <memory per CPU> : required; format: 100M, 5G etc."
    echo "       -c <CPUs per task> : optional, default: 1"
    echo "       -j <job name> : optional, default: \"noname\""
    echo "       -t <upper running time limit> : optional, default: 2h; format: hh:mm:ss"
    echo "       -r <max running jobs> : optional, default: all (if specified, should not exceed ${max_array_size})"
}

username=$(id|cut -d ' ' -f 1|sed 's/[(]/ /'|cut -d ' ' -f 2 |sed 's/[)]//')
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
runner="${script_dir}/runArrayHelper.sh"
datestring=$(date "+%Y_%b_%d-%H:%M:%S")

OPTIND=1
argfile=""
mem=""
maxtime="NA"
max_running="NA"
jobname="noname"
cpus=1
while getopts ":i:r:c:j:t:m:h" opt; do
    case $opt in
	i)
	    argfile=${OPTARG}
	    ;;
	r)
	    max_running=${OPTARG}
	    ;;
	c)
	    cpus=${OPTARG}
	    ;;
	j)
	    jobname=${OPTARG}
	    ;;
	t)
	    maxtime=${OPTARG}
	    ;;
	m)
	    mem=${OPTARG}
	    ;;
	h)
	    usage
	    exit 1
	    ;;
	\?)
	    usage
	    exit 1
	    ;;
	*)
	    usage
	    exit 1
	    ;;
    esac
done

total=$(cat $argfile|wc -l)

if [[ ${max_running} == "NA" && ${total} > ${max_array_size} ]];then
    echo "ERROR: argument file (${argfile}) is too large: (${total} lines)"
    echo "ERROR: no <max running jobs> is specified and on this cluster MaxArraySize=${max_array_size}"
    echo "ERROR: either use -r <max running jobs> with the value less than ${max_array_size} or split your argument file into chunks with #lines < ${max_array_size} each"
    exit 1
fi

if [[ ${max_running} > ${max_array_size} ]];then
    echo "ERROR: <max running jobs> is too large ( > ${max_array_size})"
    exit 1
fi

logfile=${username}_${jobname}_${datestring}.log

if [ -z ${argfile} ]; then
    echo "ERROR: no argument file name specified"
    usage
    exit 1
fi

if [ ! -f ${argfile} ]; then
    echo "ERROR: ${argfile} does not exist"
    usage
    exit 1
fi

if [ -z ${mem} ]; then
    echo "ERROR: no memory per CPU specified"
    usage
    exit 1
fi

date > ${logfile}
echo  >> ${logfile}
echo "Parameters:" >> ${logfile}
echo "commands          : $argfile" >> ${logfile}
echo "max running jobs   : $max_running" >> ${logfile}
echo "CPUs               : $cpus" >> ${logfile}
echo "job name           : $jobname" >> ${logfile}
echo "max running time   : $maxtime" >> ${logfile}
echo "mem per CPU        : $mem" >> ${logfile}
echo >> ${logfile}

CMD="sbatch --array=1-${total}"

if [[ ${max_running} != "NA" ]];then
    CMD=${CMD}"%${max_running}"
fi

if [[ ${maxtime} != "NA" ]];then
    CMD=${CMD}" --time=${maxtime}"
fi

CMD=${CMD}" --job-name=${jobname} --cpus-per-task=${cpus} ${runner} ${argfile}"
${CMD} 1 >> ${logfile} 2 >> ${logfile}
