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
    /** VALUES ARE ASSIGNED MANUALLY FOR TESTING REASONS.
    * THIS IS TEMPORARY!
    */
    int n = 840, m = 840 , mits = 50;
    double alpha = 0.8, tol = 1e-13, relax = 1.0;
    double maxAcceptableError;
    /** DONT FORGET THIS */
    double error = HUGE_VAL;
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
    clock_t start = clock(), diff;

    // Variable declarations moved outside the loop
    //macros are defined to translate the 2-D index (XX, YY) to the 1-dimensional array:
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]
    double updateVal;
    double f;
    // Coefficients
    double cx = 1.0/(deltaX*deltaX);
    double cy = 1.0/(deltaY*deltaY);
    double cc = -2.0*cx-2.0*cy-alpha;

    int x, y;
    double fX[maxXCount], fY[maxYCount];
    // Create 2 arrays with fX and fY values. These values are fixed for each iteration.
    // We are doing this in the same array because always x=y.
    for (y = 1; y < (maxYCount-1); y++)
    {
        fY[y-1] = yStart + (y-1)*deltaY;
        fX[y-1] = xStart + (y-1)*deltaX;
    }

    MPI_Comm comm;
    int size, rank, world_size, world_rank, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    //The array containing the number of processes to assign to each dimension.
    int dims[2] = {0,0};
    //logical array of size ndims specifying whether the grid is periodic ( true) or not ( false) in each dimension
    int period[] = {0,0};
    // ranking may be reordered (true) or not (false)
    int reorder = 1;
    int coords[2];

    /** MPI INIT */
    MPI_Init(NULL,NULL);
    t1 = MPI_Wtime();

    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get the rank of process
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //This routine decomposes a given number of processes
    //over a cartesian grid made of the number of dimensions specified.
    MPI_Dims_create(world_size, 2, dims);

    MPI_Get_processor_name(processor_name, &name_len);
    printf("Hello from processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, reorder, &comm);

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    //Get my coordinated in the new communicator
    MPI_Cart_coords(comm, rank, 2, coords);
    printf("[MPI process %d] I am located at (%d, %d).\n", rank, coords[0],coords[1]);

    /* Iterate as long as it takes to meet the convergence criterion */
    while (iterationCount < maxIterationCount && error > maxAcceptableError)
    {
        //printf("Iteration %i\n", iterationCount);
        error = 0.0;
        // Start of iterations
        for (y = 1; y < (maxYCount-1); y++)
        {
            for (x = 1; x < (maxXCount-1); x++)
            {
                updateVal = (	(SRC(x-1,y) + SRC(x+1,y))*cx +
                                 (SRC(x,y-1) + SRC(x,y+1))*cy +
                                 SRC(x,y)*cc
                                 -(-alpha*(1.0-fX[x]*fX[x])*(1.0-fY[y]*fY[y]) - 2.0*(1.0-fX[x]*fX[x]) - 2.0*(1.0-fY[y]*fY[y])))
                                         /cc;
                DST(x,y) = SRC(x,y) - omega*updateVal;
                error += updateVal*updateVal;
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
    printf( "Rank %d: Iterations=%3d Elapsed MPI Wall time is %f\n", rank, iterationCount, t2 - t1 );
    MPI_Finalize();


    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Rank %d: Time taken %d seconds %d milliseconds\n", rank, msec/1000, msec%1000);

    // u_old holds the solution after the most recent buffers swap
    double absoluteError = checkSolution(xLeft, yBottom,
                                         n+2, m+2,
                                         u_old,
                                         deltaX, deltaY,
                                         alpha);
    if (rank == 0) {
        printf("Residual %g\n", error);
        printf("The error of the iterative solution is %g\n", absoluteError);
    }


    return 0;
}
