#include "../include/common/read_input.h"
#include <stdio.h>

/**
 * Disable the following compiler warning for this file only:
 *  ignoring return value of ‘scanf’, declared with attribute warn_unused_result [-Wunused-result]
 */
#pragma GCC diagnostic ignored "-Wunused-result"

void read_input(
    int *n, int *m,double *alpha, double *relax, double *tol, int *mits, bool also_print_it)
{
    // Input n,m - grid dimension in x,y direction.
    scanf("%d,%d", n, m);

    // Input alpha - Helmholtz constant.
    scanf("%lf", alpha);

    // Input relax - successive over-relaxation parameter.
    scanf("%lf", relax);

    // Input tol - error tolerance for the iterative solver.
    scanf("%lf", tol);

    // Input mits - maximum solver iterations.
    scanf("%d", mits);

    if (also_print_it == true)
    {
        printf("-> %d, %d, %g, %g, %g, %d\n", *n, *m, *alpha, *relax, *tol, *mits);
    }
}
