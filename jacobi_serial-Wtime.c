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
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

/*************************************************************
 * Performs one iteration of the Jacobi method and computes
 * the residual value.
 *
 * NOTE: u(0,*), u(maxXCount-1,*), u(*,0) and u(*,maxYCount-1)
 * are BOUNDARIES and therefore not part of the solution.
 
 moved to body of for
 double one_jacobi_iteration(double xStart, double yStart,
                            int maxXCount, int maxYCount,
                            double *src, double *dst,
                            double deltaX, double deltaY,
                            double alpha, double omega)
*************************************************************/

/**********************************************************
 * Checks the error between numerical and exact solutions
 **********************************************************/
double checkSolution(double xStart, double yStart,
                     int maxXCount, int maxYCount,
                     double *u,
                     double deltaX, double deltaY,
                     double alpha)
{
#define U(XX,YY) u[(YY)*maxXCount+(XX)]
    int x, y;
    double fX, fY;
    double localError, error = 0.0;

    for (y = 1; y < (maxYCount-1); y++)
    {
        fY = yStart + (y-1)*deltaY;
        for (x = 1; x < (maxXCount-1); x++)
        {
            fX = xStart + (x-1)*deltaX;
            localError = U(x,y) - (1.0-fX*fX)*(1.0-fY*fY);
            error += localError*localError;
        }
    }
    return sqrt(error)/((maxXCount-2)*(maxYCount-2));
}


int main(int argc, char **argv)
{
    int n = 840, m = 840 , mits = 50;
    double alpha = 0.8, tol = 1e-13, relax = 1.0;
    double maxAcceptableError;
    double error;
    double *u, *u_old, *tmp;
    int allocCount;
    int iterationCount, maxIterationCount;
    double t1, t2;

//    printf("Input n,m - grid dimension in x,y direction:\n");
//    scanf("%d,%d", &n, &m);
//    printf("Input alpha - Helmholtz constant:\n");
//    scanf("%lf", &alpha);
//    printf("Input relax - successive over-relaxation parameter:\n");
//    scanf("%lf", &relax);
//    printf("Input tol - error tolerance for the iterative solver:\n");
//    scanf("%lf", &tol);
//    printf("Input mits - maximum solver iterations:\n");
//    scanf("%d", &mits);


//    printf("-> %d, %d, %g, %g, %g, %d\n", n, m, alpha, relax, tol, mits);

    allocCount = (n+2)*(m+2);
    // Those two calls also zero the boundary elements
    u = 	(double*)calloc(allocCount, sizeof(double)); //reverse order
    u_old = (double*)calloc(allocCount, sizeof(double));
    
//    printf("allocCount=%d u=%p u_old=%p\n", allocCount, u, u_old);
    
    if (u == NULL || u_old == NULL)
    {
        printf("Not enough memory for two %ix%i matrices\n", n+2, m+2);
        exit(1);
    }
    maxIterationCount = mits;
    maxAcceptableError = tol;

    // Solve in [-1, 1] x [-1, 1]
    double xLeft = -1.0, xRight = 1.0;
    double yBottom = -1.0, yUp = 1.0;

    double deltaX = (xRight-xLeft)/(n-1);
    double deltaY = (yUp-yBottom)/(m-1);
    
//  params 8
    double xStart = xLeft;
    double yStart = yBottom;
    int maxXCount = n+2;
	int maxYCount = m+2; 
    double *src   = u_old;
    double *dst  = u;
    double omega  = relax;
    alpha         = alpha;
    
    iterationCount = 0;
    error = HUGE_VAL;
    clock_t start = clock(), diff;
    
    MPI_Init(NULL,NULL);
    t1 = MPI_Wtime();

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    //Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    printf(".Hello world from processor %s, rank %d out of %d processors\n", processor_name, world_rank, world_size);
    /* Iterate as long as it takes to meet the convergence criterion */
    while (iterationCount < maxIterationCount && error > maxAcceptableError)
    {    	
//        printf("Iteration %i", iterationCount);


#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]
    int x, y;
    double fX, fY;
    double error = 0.0;
    double updateVal;
    double f;
    // Coefficients
    double cx = 1.0/(deltaX*deltaX);
    double cy = 1.0/(deltaY*deltaY);
    double cc = -2.0*cx-2.0*cy-alpha;

    for (y = 1; y < (maxYCount-1); y++)
    {

        for (x = 1; x < (maxXCount-1); x++)
        {
            /***
             * Step A.2.b from 'instructions for effective MPI'
             * Remove variable assignments (fX, fY). Replace fx, fy in f with their values.
             * Did the same below with 'updateVal'.
             */

            f = -alpha*(1.0-(xStart + (x-1)*deltaX)*(xStart + (x-1)*deltaX))
                      *(1.0-(yStart + (y-1)*deltaY)*(yStart + (y-1)*deltaY))
                      - 2.0*(1.0-(xStart + (x-1)*deltaX)*(xStart + (x-1)*deltaX))
                      - 2.0*(1.0-(yStart + (y-1)*deltaY)*(yStart + (y-1)*deltaY));

            DST(x,y) = SRC(x,y) - omega*updateVal;
            error += ((SRC(x-1,y) + SRC(x+1,y))*cx +
                          (SRC(x,y-1) + SRC(x,y+1))*cy +
                          SRC(x,y)*cc - f
                     )/cc
                             *
                    ((SRC(x-1,y) + SRC(x+1,y))*cx +
                         (SRC(x,y-1) + SRC(x,y+1))*cy +
                         SRC(x,y)*cc - f
                    )/cc;
        }
    }
    error= sqrt(error)/((maxXCount-2)*(maxYCount-2));
                            
//        printf("\tError %g\n", error);
        iterationCount++;
        // Swap the buffers
        tmp = u_old;
        u_old = u;
        u = tmp;
    }

    t2 = MPI_Wtime();
    printf( "Iterations=%3d Elapsed MPI Wall time is %f\n", iterationCount, t2 - t1 ); 
    MPI_Finalize();


    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("Residual %g\n",error);

    // u_old holds the solution after the most recent buffers swap
    double absoluteError = checkSolution(xLeft, yBottom,
                                         n+2, m+2,
                                         u_old,
                                         deltaX, deltaY,
                                         alpha);
    printf("The error of the iterative solution is %g\n", absoluteError);

    return 0;
}
