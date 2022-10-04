#!/bin/bash

absolute_path=`cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && (pwd -W 2> /dev/null || pwd)`
common_scripts_path="$absolute_path/../../common/scripts"

. "$common_scripts_path/print.sh"

################################################################################################

report_output_file="./argo/reports/mpi_performance_`date '+%Y_%m_%d_%H_%M_%S'`.csv"

# Add CSV header row.
echo `$common_scripts_path/csv_create_header_row.sh`",out_size,out_np" > "$report_output_file"

# First test it for small sizes. If everything goes fine, then try the larger sizes.
#for size in 840 1680 3360 6270 13440 26880; do
for size in 840 1680 3360; do
for np in 1 4 9 16 25 36 49 64 80; do
for i in {1..3}; do

    # Create new qsub job.
    print_green "make --no-print-directory x size=$size np=$np"
    make --no-print-directory x "size=$size" "np=$np"
    job_id=`cat .latest_qsub_job_id`

    # Wait for the qsub job to finish.
    $common_scripts_path/wait_qsub_job.sh "$job_id"
    [[ $? != 0 ]] && exit 1;

    # Add CSV data row.
    job_output_file="./argo/outputs/$job_id.OU"
    csv_data_row=`$common_scripts_path/csv_create_data_row.sh < "$job_output_file"`",$size,$np"
    print_green "$csv_data_row"
    echo "$csv_data_row" >> "$report_output_file"

    # Do not spawn next job immediately.
    sleep 2
done
done
done


