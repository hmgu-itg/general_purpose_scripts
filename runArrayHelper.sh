#!/bin/bash

argfile=$1

echo "RUNNING ON" $(hostname)
argline=$(head -n $SLURM_ARRAY_TASK_ID $argfile|tail -n 1)
read -r -a arguments <<< "${argline}"
#singularity exec -B /compute/Genomics /compute/Genomics/containers/worker_1.8s "${arguments[@]}"
"${arguments[@]}"

exit 0
