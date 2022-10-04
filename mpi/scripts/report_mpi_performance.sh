#!/bin/bash

absolute_path=`cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && (pwd -W 2> /dev/null || pwd)`
common_scripts_path="$absolute_path/../../common/scripts"

. "$common_scripts_path/print.sh"

################################################################################################

report_output_file="./argo/reports/mpi_performance_`date '+%Y_%m_%d_%H_%M_%S'`.csv"

# Add CSV header row.
echo `$common_scripts_path/csv_create_header_row.sh`",out_np" > "$report_output_file"

for size in 840 1680 3360 6270 13440 26880; do
for np in 1 4 9 16 25 36 49 64 80; do
  if [[ $size == 13440 ]]; then
      if [[ $np < 9 ]]; then
          repeats=2
      fi
  elif [[ $size == 26880 ]]; then
      if [[ $np < 9 ]]; then
          repeats=1
      elif [[ $np < 25 ]]; then
          repeats=2
      else
          repeats=3
      fi
  else
      repeats=3
  fi
  for i in `seq 1 $repeats`; do

    # Create new qsub job.
    print_green "make --no-print-directory x size=$size np=$np"
    make --no-print-directory x "size=$size" "np=$np"
    job_id=`cat .latest_qsub_job_id`

    # Wait for the qsub job to finish.
    $common_scripts_path/wait_qsub_job.sh "$job_id"
    [[ $? != 0 ]] && exit 1;

    # Add CSV data row.
    job_output_file="./argo/outputs/$job_id.OU"
    csv_data_row=`$common_scripts_path/csv_create_data_row.sh < "$job_output_file"`",$np"
    print_green "$csv_data_row"
    echo "$csv_data_row" >> "$report_output_file"

    # Do not spawn next job immediately.
    sleep 2
done
done
done


