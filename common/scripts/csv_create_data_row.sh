#!/bin/bash

read_input_line() {
    IFS= read -r line
    echo $line
}

# Parses a job's output (which is supposed to be the same for all modules: sequential, mpi, openmp, etc)
# and reports the information of interest, comma separated.
#
# Example job output and input to this script:
#<- 840, 840, 0.8, 1, 1e-13, 50
#-> Iterations: 50
#-> MPI_Wtime: 1.229024 secs, clock: 1.222 secs
#-> Residual: 5.40181e-09
#-> Iterative solution error: 0.000633904
# 
# Example script's output based on the above input:
# 840,840,0.8,1,1e-13,50,50,0.799028,0.797,5.40181e-09,0.000633904,1,0

read_input_line | sed -e 's/,//g' | awk 'BEGIN{ OFS=","; ORS="," } {print $2,$3,$4,$5,$6,$7}'
read_input_line | awk 'BEGIN{ ORS="," } {print $3}'
read_input_line | awk 'BEGIN{ OFS=","; ORS="," } {print $3,$6}'
read_input_line | awk 'BEGIN{ ORS="," } {print $3}'
read_input_line | awk 'BEGIN{ ORS="" } {print $5}'


