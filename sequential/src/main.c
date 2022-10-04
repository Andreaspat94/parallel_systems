/************************************************************
 * Program to solve a finite difference
 * discretization of the screened Poisson equation:
 * (d2/dx2)u + (d2/dy2)u - alpha u = f
 * with zero Dirichlet boundary condition using the iterative
 * Jacobi method with overrelaxation.
 *
 * RHS (source) function
 *   f(x,y) = -alpha*(1-x^2)(1-y^2)-2*[(1-x^2)+(1-y^2)]
 *
 * Analytical solution to the PDE
 *   u(x,y) = (1-x^2)(1-y^2)
 *
 * Current Version: Christian Iwainsky, RWTH Aachen University
 * MPI C Version: Christian Terboven, RWTH Aachen University, 2006
 * MPI Fortran Version: Dieter an Mey, RWTH Aachen University, 1999 - 2005
 * Modified: Sanjiv Shah,        Kuck and Associates, Inc. (KAI), 1998
 * Author:   Joseph Robicheaux,  Kuck and Associates, Inc. (KAI), 1998
 *
 * Unless READ_INPUT is defined, a meaningful input dataset is used (CT).
 *
 * Input : n     - grid dimension in x direction
 *         m     - grid dimension in y direction
 *         alpha - constant (always greater than 0.0)
 *         tol   - error tolerance for the iterative solver
 *         relax - Successice Overrelaxation parameter
 *         mits  - maximum iterations for the iterative solver
 *
 * On output
 *       : u(n,m)       - Dependent variable (solution)
 *       : f(n,m,alpha) - Right hand side function
 *
 *************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
#include "common/input.h"
#include "common/prints.h"
#include "common/timing.h"
#include "common/allocate_grid.h"
#include "common/check_solution.h"
#include "precalculations_t.h"
#include "jacobi_iteration_params_t.h"
#include "jacobi_iteration_original.h"
#include "jacobi_iteration_opt1.h"
#include "jacobi_iteration_opt1x.h"
#include "jacobi_iteration_opt2.h"
#include "jacobi_iteration_opt2x.h"
#include "jacobi_iteration_opt3.h"
#include "jacobi_iteration_opt3x.h"

int main(int argc, char **argv)
{
    int rank;
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ////////////////////////////////////////////////////////////////////////////////////////////////

    input_t input;
    input_read_parallel(&input, rank, MPI_COMM_WORLD);

    if (rank == 0)
        print_input(&input);

    int n = input.n;
    int m = input.m;
    double alpha = input.alpha;
    double relax = input.relax;
    double max_acceptable_error = input.max_acceptable_error;
    int max_iteration_count = input.max_iteration_count;

    // Solve in [-1, 1] x [-1, 1]
    double xLeft = -1.0, xRight = 1.0;
    double yBottom = -1.0, yUp = 1.0;

    double deltaX = (xRight-xLeft)/(n-1);
    double deltaY = (yUp-yBottom)/(m-1);

    ////////////////////////////////////////////////////////////////////////////////////////////////

    void (*jacobi_precalculate)(double, double, int, int, double, double, double,
        precalculations_t*);
    double (*jacobi_iteration)(jacobi_iteration_params_t*, precalculations_t*);
    {
        if (argc > 2)
        {
            fprintf(stderr, "No more than 1 argument.");
            exit(1);
        }
        else if (argc == 2)
        {
            const char *arg = argv[1];
            if (!strcmp("-O0", arg))
            {
                jacobi_precalculate = &jacobi_precalculate_original;
                jacobi_iteration    = &jacobi_iteration_original;
            }
            else if (!strcmp("-O1", arg))
            {
                jacobi_precalculate = &jacobi_precalculate_opt1;
                jacobi_iteration    = &jacobi_iteration_opt1;
            }
            else if (!strcmp("-O1x", arg))
            {
                jacobi_precalculate = &jacobi_precalculate_opt1x;
                jacobi_iteration    = &jacobi_iteration_opt1x;
            }
            else if (!strcmp("-O2", arg))
            {
                jacobi_precalculate = &jacobi_precalculate_opt2;
                jacobi_iteration    = &jacobi_iteration_opt2;
            }
            else if (!strcmp("-O2x", arg))
            {
                jacobi_precalculate = &jacobi_precalculate_opt2x;
                jacobi_iteration    = &jacobi_iteration_opt2x;
            }
            else if (!strcmp("-O3", arg))
            {
                jacobi_precalculate = &jacobi_precalculate_opt3;
                jacobi_iteration    = &jacobi_iteration_opt3;
            }
            else if (!strcmp("-O3x", arg))
            {
                jacobi_precalculate = &jacobi_precalculate_opt3x;
                jacobi_iteration    = &jacobi_iteration_opt3x;
            }
            else
            {
                fprintf(stderr, "Valid arguments: \"-O0\" (default), \"-O1\", \"-O1x\", \"-O2\", \"-O2x\", \"-O3\", \"-O3x\"");
                exit(1);
            }
        }
        else
        {
            jacobi_precalculate = &jacobi_precalculate_original;
            jacobi_iteration    = &jacobi_iteration_original;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////

    double *u, *u_old;
    allocate_grid(n+2, m+2, &u, &u_old);

    precalculations_t precalculations;
    (*jacobi_precalculate)(xLeft, yBottom, n, m, deltaX, deltaY, alpha, &precalculations);

    jacobi_iteration_params_t jacobi_iteration_params = {
        xLeft, yBottom, n+2, m+2, u_old, u, deltaX, deltaY, alpha, relax
    };

    int iteration_count = 0;
    double error = HUGE_VAL;

    times_t times;
    times_begin(&times);

    // Iterate as long as it takes to meet the convergence criterion.
    while (iteration_count < max_iteration_count && error > max_acceptable_error)
    {
        // printf("Iteration %i", iterationCount);

        error = jacobi_iteration(&jacobi_iteration_params, &precalculations);

        // printf("\tError %g\n", error);

        swap_buffers(&jacobi_iteration_params);

        iteration_count++;
    }

    times_end(&times);
    times_reduce_max(&times, 0, MPI_COMM_WORLD);

    free(precalculations.f);
    precalculations.f = NULL;

    // "jacobi_iteration_params.src" holds the solution after the most recent buffers swap.
    double absolute_error = check_solution(
        xLeft, yBottom,
        n+2, m+2,
        jacobi_iteration_params.src,
        deltaX, deltaY);
    absolute_error = sqrt(absolute_error)/(n*m);

    if (rank == 0)
        print_output(iteration_count, &times, error, absolute_error);

    MPI_Finalize();

    return 0;
}
