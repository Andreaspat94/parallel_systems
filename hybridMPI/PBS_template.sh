#!/bin/bash

# JobName #
#PBS -N openMPI

#Which Queue to use #
#PBS -q N10C80

# Max Wall time, Example 1 Minute #
#PBS -l walltime=00:20:00

# How many nodes and tasks per node
#PBS -l select=%SELECT%:ncpus=8:mpiprocs=%MPIPROCS%:ompthreads=%OMPTHREADS%:mem=15G

#Change Working directory to SUBMIT directory
cd $PBS_O_WORKDIR

# Run executable #
mpirun -np %NP% ./argo/bin/%TARGET%.x < ./argo/inputs/input_%SIZE%.txt


