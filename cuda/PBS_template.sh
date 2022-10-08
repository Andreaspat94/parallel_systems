#!/bin/bash

# JobName
#PBS -N cuda

# Which Queue to use
#PBS -q GPUq

# Max Wall time, Example 1 Minute
#PBS -l walltime=00:01:00

# How many nodes, cpus/node, mpiprocs/node and threads/mpiprocess
#PBS -l select=1:ncpus=2:ompthreads=2:ngpus=%NGPUS% -lplace=excl

# Change Working directory to SUBMIT directory
cd $PBS_O_WORKDIR

# Run executable.
./argo/bin/jacobi_gpu.x < ./argo/inputs/input_%SIZE%.txt

