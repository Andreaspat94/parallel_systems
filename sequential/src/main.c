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
#include "common/read_input.h"
#include "common/allocate_grid.h"
#include "common/check_solution.h"
#include "jacobi_iteration_original.h"
#include "jacobi_iteration_opt1.h"
#include "jacobi_iteration_opt2.h"

int main(int argc, char **argv)
{
    int n, m;
    double alpha, relax;
    double maxAcceptableError;
    int maxIterationCount;
    read_input(&n, &m, &alpha, &relax, &maxAcceptableError, &maxIterationCount, true);

    // Solve in [-1, 1] x [-1, 1]
    double xLeft = -1.0, xRight = 1.0;
    double yBottom = -1.0, yUp = 1.0;

    double deltaX = (xRight-xLeft)/(n-1);
    double deltaY = (yUp-yBottom)/(m-1);

    ////////////////////////////////////////////////////////////////////////////////////////////////

    void (*jacobi_precalculate)(double, double, int, int, double, double, double,
        precalculations_t*);
    double (*jacobi_iteration)(double, double, int, int, const double*, double*, double, double,
        double, double, precalculations_t*);
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
            else if (!strcmp("-O2", arg))
            {
                jacobi_precalculate = &jacobi_precalculate_opt2;
                jacobi_iteration    = &jacobi_iteration_opt2;
            }
            else
            {
                fprintf(stderr, "Valid arguments: \"-O0\" (default), \"-O1\"");
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

    precalculations_t precalculations;
    (*jacobi_precalculate)(xLeft, yBottom, n, m, deltaX, deltaY, alpha, &precalculations);

    double *u, *u_old;
    allocate_grid(n, m, &u, &u_old);

    int iterationCount = 0;
    double error = HUGE_VAL;

    clock_t clock1 = clock();
    MPI_Init(NULL,NULL);
    double t1 = MPI_Wtime();

    // Iterate as long as it takes to meet the convergence criterion.
    while (iterationCount < maxIterationCount && error > maxAcceptableError)
    {
        // printf("Iteration %i", iterationCount);

        error = jacobi_iteration(
            xLeft, yBottom,
            n+2, m+2,
            u_old, u,
            deltaX, deltaY,
            alpha, relax,
            &precalculations
        );

        // printf("\tError %g\n", error);

        // Swap the buffers
        double *tmp = u_old;
        u_old = u;
        u = tmp;

        iterationCount++;
    }

    double t2 = MPI_Wtime();
    MPI_Finalize();
    clock_t clock2 = clock();
    clock_t msec = (clock2 - clock1) * 1000 / CLOCKS_PER_SEC;
    printf("Iterations=%3d Elapsed MPI Wall time is %f\n", iterationCount, t2 - t1);
    printf("Time taken %ld seconds %ld milliseconds\n", msec/1000, msec%1000);

    printf("Residual %g\n",error);

    // u_old holds the solution after the most recent buffers swap.
    double absoluteError = check_solution(
        xLeft, yBottom,
        n+2, m+2,
        u_old,
        deltaX, deltaY);

    printf("The error of the iterative solution is %g\n", absoluteError);

    // TODO: free "precalculations.f".

    return 0;
}
