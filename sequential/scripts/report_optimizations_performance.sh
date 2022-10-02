#!/bin/bash

absolute_path=`cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && (pwd -W 2> /dev/null || pwd)`
common_scripts_path="$absolute_path/../../common/scripts"

. "$common_scripts_path/print.sh"

################################################################################################

report_output_file="./argo/reports/optimizations_performance_`date '+%Y_%m_%d_%H_%M_%S'`.csv"

# Add CSV header row.
echo `$common_scripts_path/csv_create_header_row.sh`",out_np,out_opt" > "$report_output_file"

for np in {1..8}; do
for opt in 0 1 1x 2 2x 3 3x; do

    # Create new qsub job.
    print_green "make --no-print-directory x np=$np opt=$opt"
    make --no-print-directory x "np=$np" "opt=$opt"
    job_id=`cat .latest_qsub_job_id`

    # Wait for the qsub job to finish.
    $common_scripts_path/wait_qsub_job.sh "$job_id"
    [[ $? != 0 ]] && exit 1;

    # Add CSV data row.
    job_output_file="./argo/outputs/$job_id.OU"
    csv_data_row=`$common_scripts_path/csv_create_data_row.sh < "$job_output_file"`",$np,$opt"
    print_green "$csv_data_row"
    echo "$csv_data_row" >> "$report_output_file"

    # Do not spawn next job immediately.
    sleep 2
done
done


