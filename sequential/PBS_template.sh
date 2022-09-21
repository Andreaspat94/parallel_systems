#!/bin/bash

# JobName #
#PBS -N sequential

#Which Queue to use #
#PBS -q N10C80

# Max Wall time, Example 1 Minute #
#PBS -l walltime=00:20:00

# How many nodes and tasks per node
#PBS -l select=1:ncpus=8:mem=2G

#Change Working directory to SUBMIT directory
cd $PBS_O_WORKDIR

# Run executable #
mpirun -np 1 ./argo/bin/%TARGET%.x < ./argo/inputs/input_%SIZE%.txt


