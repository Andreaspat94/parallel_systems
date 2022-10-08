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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "common/check_solution.h"
#include "jacobi_gpu.h"

int main(int argc, char **argv)
{
    int n, m;
    double alpha, relax;
    double max_acceptable_error;
    int max_iteration_count;
    
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
    // Input n,m - grid dimension in x,y direction.
    scanf("%d,%d", &n, &m);
    // Input alpha - Helmholtz constant.
    scanf("%lf", &alpha);
    // Input relax - successive over-relaxation parameter.
    scanf("%lf", &relax);
    // Input tol - error tolerance for the iterative solver.
    scanf("%lf", &max_acceptable_error);
    // Input mits - maximum solver iterations.
    scanf("%d", &max_iteration_count);
#pragma GCC diagnostic pop
    
    printf("<- %d, %d, %g, %g, %g, %d\n",
           n, m, alpha, relax, max_acceptable_error, max_iteration_count);

    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Solve in [-1, 1] x [-1, 1]
    double xLeft = -1.0, xRight = 1.0;
    double yBottom = -1.0, yUp = 1.0;

    double deltaX = (xRight-xLeft)/(n-1);
    double deltaY = (yUp-yBottom)/(m-1);

    int maxXCount = n*2;
    int maxYCount = m*2;

    ////////////////////////////////////////////////////////////////////////////////////////////////

    auto *u     = (double*)calloc(maxXCount*maxYCount, sizeof(double));
    auto *u_old = (double*)calloc(maxXCount*maxYCount, sizeof(double));
    if (u == nullptr || u_old == nullptr)
    {
        printf("Not enough memory for two %ix%i matrices\n", maxXCount, maxYCount);
        exit(1);
    }

    double error;
    int iteration_count;
    float elapsedTime;

    jacobi_gpu(u_old, u,
               n+2, m+2,
               xLeft, yBottom,
               deltaX, deltaY,
               alpha, relax,
               max_iteration_count, max_acceptable_error,
               &iteration_count, &error, &elapsedTime);

    // "u_old" holds the solution after the most recent buffers swap.
    double absolute_error = check_solution(
        xLeft, yBottom,
        n+2, m+2,
        u_old,
        deltaX, deltaY);
    absolute_error = sqrt(absolute_error)/(n*m);

    printf("-> Iterations: %d\n", iteration_count);
    printf("-> Time: %d.%03d secs\n", (int)elapsedTime/1000, ((int)elapsedTime)%1000);
    printf("-> Residual: %g\n", error);
    printf("-> Iterative solution error: %g\n", absolute_error);

    free(u);
    free(u_old);

    return 0;
}
