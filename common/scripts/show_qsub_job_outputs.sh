#!/bin/bash

fullpath=`cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && (pwd -W 2> /dev/null || pwd)`

. "$fullpath/print.sh"

################################################################################################

if [[ $# -eq 0 ]]; then
    # Specific job id was not given. Trying ".latest_qsub_job_id" file.
    if [ ! -f .latest_qsub_job_id ]; then
        print_error 'No ".latest_qsub_job_id" file found. Aborting...'
        exit 1
    fi
    job_id="`cat .latest_qsub_job_id`"
elif [[ "$1" =~ ^[0-9]+\.argo$ ]]; then
    # Specific job id given.
    job_id="$1"
elif [[ "$1" =~ ^[0-9]+$ ]]; then
    # Specific job id given, but only the number part.
    job_id="$1.argo"
else
    print_error 'Given argument must be a qsub job id, either its number part, or the whole id. Aborting...'
    exit 1
fi

################################################################################################

function print_header {
    local type="$1"
    local file="$2"
    local file_creation_date="`stat \"$file\" | tail -n 1 | sed -e 's/Birth://' -e 's/\.[0-9]\{9\}//g' | xargs`"    

    python -c 'print("~" * 80)'
    python -c "print('~' * 4 + ' ($type | $file | $file_creation_date)')"
    python -c 'print("~" * 80)'
}

function print_body {
    local type="$1"
    local file="$2"

    if [[ "$type" == "stdout" ]]; then
        print_green "$(cat $file)"
    elif [[ "$type" == "stderr" ]]; then
        print_red "`cat $file`"
    else
        cat $file
    fi
}

function print_block {
    local type="$1"
    local file="$2"

    if [ ! -f "$file" ]; then
        # File does not exist.
        return 1
    elif [[ "`wc -c \"$file\" | cut -d \" \" -f1`" == "0" ]]; then
        # File is empty.
        return 0
    fi

    echo
    print_header $type $file
    print_body $type $file
}

print_block stdout "./argo/outputs/$job_id.OU"
if [ "$?" == "1" ]; then
    echo
    print_error "Did not find any files for job \"$job_id\" in \"./argo/outputs/\""
    echo
    exit 1
fi
print_block stderr "./argo/outputs/$job_id.ER"
echo


