#ifndef SEQUENTIAL_TIMING_H
#define SEQUENTIAL_TIMING_H

#include <stdio.h>
#include <mpi.h>
#include <time.h>

typedef struct {
    double mpi_wall[2];
    clock_t clock[2];
} times_t;

void times_begin(times_t *times);
void times_end(times_t *times);
double times_mpi_wtime_dt(const times_t *times);
clock_t times_clock_dt(const times_t *times);
void times_reduce_max(times_t *times, int rank_dst, MPI_Comm comm);

#endif //SEQUENTIAL_TIMING_H
