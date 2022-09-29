#ifndef SEQUENTIAL_INPUT_H
#define SEQUENTIAL_INPUT_H

#include <mpi.h>

typedef struct {
    int n;
    int m;
    double alpha;
    double relax;
    double max_acceptable_error;
    int max_iteration_count;
} input_t;

void input_read(input_t *input);
void input_read_parallel(input_t *input, int rank, MPI_Comm comm);

#endif //SEQUENTIAL_INPUT_H
