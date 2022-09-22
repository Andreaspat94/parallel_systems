#!/bin/bash

absolute_path=`cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && (pwd -W 2> /dev/null || pwd)`

. "$absolute_path/print.sh"

################################################################################################

if [[ ! -f .latest_qsub_job_id ]]; then
    print_error 'Could not rename PBS configuration file. No ".latest_qsub_job_id" file found.'
    exit 1
fi

job_id=`cat .latest_qsub_job_id`

if [[ ! $job_id =~ ^[0-9]+\.argo$ ]]; then
    print_error 'Could not rename PBS configuration file. Invalid job id format inside ".latest_qsub_job_id".'
    exit 1
fi

file_to_rename="$1"
mv "$file_to_rename" "./argo/inputs/$job_id.pbs_config"


