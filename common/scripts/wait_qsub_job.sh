#!/bin/bash

absolute_path=`cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && (pwd -W 2> /dev/null || pwd)`

. "$absolute_path/print.sh"

################################################################################################

job_id=$1

if [[ ! $1 =~ ^[0-9]+$ ]] && [[ ! $1 =~ ^[0-9]+\.argo$ ]]; then
    print_error 'Given argument must be a qsub job id, either its number part, or the whole id. Aborting...'
    exit 1
fi

non_running_efforts=0

while true; do
    line=`qstat -H "$job_id" | tail -1`
    state=`echo $line | awk '{print $10}'`
    print_yellow $line

    if [[ $state == F ]]; then
        break
    elif [[ $state == R ]]; then
        sleep 2
    elif [[ $non_running_efforts == 3 ]]; then
        print_error "Job '$job_id' stuck in queue for too long. Aborting..."
        exit 1
    else
        non_running_efforts=$((non_running_efforts+1))
        sleep 2
    fi
done

# Wait a bit in case the output creation somehow delays.
sleep 2

output_file="./argo/outputs/$job_id.OU"
if [[ ! -f $output_file ]]; then
    print_error "Output file for job '$job_id' was not found. Aborting..."
    exit 1
fi


