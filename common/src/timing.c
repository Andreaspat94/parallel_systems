#include "../include/common/timing.h"

void times_begin(times_t *times)
{
    times->mpi_wall[0] = MPI_Wtime();
    times->clock[0] = clock();
}

void times_end(times_t *times)
{
    times->mpi_wall[1] = MPI_Wtime();
    times->clock[1] = clock();
}

double times_mpi_wtime_dt(const times_t *times)
{
    return times->mpi_wall[1] - times->mpi_wall[0];
}

clock_t times_clock_dt(const times_t *times)
{
    return (times->clock[1] - times->clock[0]) * 1000 / CLOCKS_PER_SEC;
}

void times_reduce_max(times_t *times, int rank_dst, MPI_Comm comm)
{
    times_t temp = *times;
    MPI_Reduce(&temp.mpi_wall, &times->mpi_wall, 1, MPI_DOUBLE, MPI_MAX, rank_dst, comm);
    MPI_Reduce(&temp.clock, &times->clock, 1, MPI_LONG, MPI_MAX, rank_dst, comm);
}
