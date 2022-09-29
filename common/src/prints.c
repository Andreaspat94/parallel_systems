#include "../include/common/prints.h"
#include <stdio.h>

void print_input(const input_t *input)
{
    printf("<- %d, %d, %g, %g, %g, %d\n",
           input->n, input->m,
           input->alpha, input->relax,
           input->max_acceptable_error, input->max_iteration_count);
    fflush(stdout);
}

void print_iterations(int iteration_count)
{
    printf("-> Iterations: %d\n", iteration_count);
}

void print_times(const times_t *times)
{
    printf("-> MPI_Wtime: %f secs, clock: %ld.%03ld secs\n",
           times_mpi_wtime_dt(times),
           times_clock_dt(times)/1000,
           times_clock_dt(times)%1000);
}

void print_residual(double error)
{
    printf("-> Residual: %g\n", error);
}

void print_iterative_solution_error(double absolute_error)
{
    printf("-> Iterative solution error: %g\n", absolute_error);
}

void print_output(int iteration_count, const times_t *times, double error, double absolute_error)
{
    print_iterations(iteration_count);
    print_times(times);
    print_residual(error);
    print_iterative_solution_error(absolute_error);
    fflush(stdout);
}
