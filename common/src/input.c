#include "../include/common/input.h"
#include <stdio.h>

void input_read(input_t *input)
{
    // Disable the following compiler warning for this function only:
    //  ignoring return value of ‘scanf’, declared with attribute warn_unused_result [-Wunused-result]
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"

    // Input n,m - grid dimension in x,y direction.
    scanf("%d,%d", &input->n, &input->m);

    // Input alpha - Helmholtz constant.
    scanf("%lf", &input->alpha);

    // Input relax - successive over-relaxation parameter.
    scanf("%lf", &input->relax);

    // Input tol - error tolerance for the iterative solver.
    scanf("%lf", &input->max_acceptable_error);

    // Input mits - maximum solver iterations.
    scanf("%d", &input->max_iteration_count);

#pragma GCC diagnostic pop
}

void input_read_parallel(input_t *input, int rank, MPI_Comm comm)
{
    if (rank == 0)
    {
        input_read(input);
    }

    MPI_Bcast(input, sizeof(input_t), MPI_BYTE, 0, comm);
}
