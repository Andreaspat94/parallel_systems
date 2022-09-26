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
#include "common/read_input.h"
#include "common/allocate_grid.h"
#include "common/check_solution.h"

/**
 * Changes:
 *
 * - Require that NP is either 80 or a perfect square.
 *
 * - Read problem parameters from input file instead of have them hardcoded.
 *      use a custom struct "bcast_input_t" just for faster sharing of these values.
 *
 * - Introduce "comm_t" struct to encapsulate the 3 related variables (i.e. MPI_Comm, size, rank).
 *
 * - MPI_Cart_shift(comm, 1, 1, &south, &north); => changed to first "north" then "south".
 *
 * - NOT ----------------------Swapped directions in MPI_Cart_shift: !!!!! NOT !!!!!
 *      MPI_Cart_shift(comm, 1, 1, &south, &north); => vertical direction is "0".
 *      MPI_Cart_shift(comm, 0, 1, &west, &east);   => horizontal direction is "1".
 */

typedef struct {
    int n;
    int m;
    double alpha;
    double relax;
    double max_acceptable_error;
    int max_iteration_count;
} bcast_input_t;

typedef struct {
    MPI_Comm id;
    int rank;
    int size;
} comm_t;

typedef struct {
    int north;
    int east;
    int south;
    int west;
} neighbor_ranks_t;

bool is_perfect_square(int number)
{
    double root = sqrt(number);
    return root == floor(root);
}

int main(int argc, char **argv)
{
    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Initialize MPI and collect MPI_COMM_WORLD-related info.

    MPI_Init(NULL,NULL);
    MPI_Barrier(MPI_COMM_WORLD); // TODO: why barrier here?

    comm_t comm_world = { MPI_COMM_WORLD };
    MPI_Comm_size(comm_world.id, &comm_world.size);
    MPI_Comm_rank(comm_world.id, &comm_world.rank);

    // Require that the total number of processes ∈ {1,4,9,16,25,36,49,64,80}.
    if (comm_world.size != 80 && !(is_perfect_square(comm_world.size) && comm_world.size < 80))
    {
//        MPI_Abort(comm_world.id, 1);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Read input values from stdin (only main process) and broadcast them to all processes.

    int n_global, m_global;
    double alpha, relax;
    double max_acceptable_error;
    int max_iteration_count;

    bcast_input_t input;

    if (comm_world.rank == 0)
    {
        read_input(&n_global, &m_global, &alpha, &relax,
                   &max_acceptable_error, &max_iteration_count, true);
    }

    MPI_Bcast(&input, sizeof(bcast_input_t), MPI_BYTE, 0, comm_world.id);

    if (comm_world.rank != 0)
    {
        n_global = input.n;
        m_global = input.m;
        alpha = input.alpha;
        relax = input.relax;
        max_acceptable_error = input.max_acceptable_error;
        max_iteration_count = input.max_iteration_count;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Create a cartesian topology with reordered process ranks and retrieve the per-process
    /// coordinates inside the 2D topology.
    ///
    /// https://www.open-mpi.org/doc/v3.1/man3/MPI_Dims_create.3.php
    /// https://www.open-mpi.org/doc/v3.1/man3/MPI_Cart_create.3.php
    /// https://www.open-mpi.org/doc/v3.0/man3/MPI_Cart_shift.3.php

    int dims[] = {0,0};     // Number of processes per dimension (assigned by "MPI_Dims_create").

    // Assign number of processes to each dimension.
    MPI_Dims_create(comm_world.size, 2, dims);

    comm_t comm_cart;       // The cartesian topology's new communicator data (later assigned).
    int periods[] = {0,0};  // Whether the grid is periodic or not (i.e. wraps around) per dimension.
    int reorder = 1;        // Whether to reorder process ranks or not.
                            // Rank reordering may improve performance in some MPI implementations.

    // Create the cartesian topology & retrieve its new communicators' related info.
    MPI_Cart_create(comm_world.id, 2, dims, periods, reorder, &comm_cart.id);
    MPI_Comm_size(comm_cart.id, &comm_cart.size);
    MPI_Comm_rank(comm_cart.id, &comm_cart.rank);

    // Get my coordinates inside the cartesian topology.
    int coords[2];
    MPI_Cart_coords(comm_cart.id, comm_cart.rank, 2, coords);

    neighbor_ranks_t ranks; // Stencil/Neighbor process ranks inside the topology.

    // direction = 0 or 1 corresponding to the two dimensions x,y.
    MPI_Cart_shift(comm_cart.id, 0, 1, &ranks.west,  &ranks.east);
    MPI_Cart_shift(comm_cart.id, 1, 1, &ranks.north, &ranks.south);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Create MPI datatype for halo points.

    int n = n_global / dims[0]; // Local row size (i.e. number of columns).
    int m = m_global / dims[1]; // Local column size (i.e. number of rows).

    MPI_Datatype row;
    MPI_Type_contiguous(n, MPI_DOUBLE, &row);
    MPI_Type_commit(&row);

    MPI_Datatype column;
    MPI_Type_vector(m, 1, n, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
//    {
//        MPI_Comm comm = comm_cart.id;
//        int rank = comm_cart.rank;
//        int size = comm_cart.size;
//        char _;
//        if (rank > 0)
//        {
//            MPI_Recv(&_, 1, MPI_BYTE, rank-1, 0, comm, MPI_STATUS_IGNORE);
//            if (rank < size-1)
//                MPI_Send(&_, 1, MPI_BYTE, rank+1, 0, comm);
//        }
//
//        printf("C[%2d/%d]: neighbor_ranks(W:%2d,E:%2d,N:%2d,S:%2d), my_coords:(%d,%d)\n",
//               rank, size,
//               ranks.west, ranks.east, ranks.north, ranks.south,
//               coords[0], coords[1]);fflush(stdout);
//
//        if (rank == 0 && size > 1)
//        {
//            printf("====> %d,%d\n", dims[0], dims[1]);fflush(stdout);
//            MPI_Send(&_, 1, MPI_BYTE, 1, 0, comm);
//        }
//
//        MPI_Comm_free(&comm_cart.id);
//        MPI_Barrier(comm_world.id);
//
//        MPI_Finalize();
//        return 0;
//    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Calculate which rectangular parts of [-1,1]x[-1,1] will be assigned to each process.

    // Calculate the part of [-1, 1] that will be assigned to this process.
    double xLeft  = -1.0,
           xRight =  1.0,
           xSlice = (xRight - xLeft) / dims[0];
    xLeft  = -1.0 + xSlice * coords[0];
    xRight = xLeft + xSlice;

    // Calculate the part of [-1, 1] that will be assigned to this process.
    double yBottom = -1.0,
           yUp     =  1.0,
           ySlice  = (yUp - yBottom) / dims[1];
    yBottom = -1.0 + ySlice * coords[1];
    yUp     = yBottom + ySlice;

    double deltaX = (xRight-xLeft)/(n-1);
    double deltaY = (yUp-yBottom)/(m-1);
    double xStart = xLeft;
    double yStart = yBottom;
    int maxXCount = n+2;
    int maxYCount = m+2;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Make any jacobi-iteration precalculations.

    double cx, cy, cc;
    double *fX, *fY;
    {
        // Coefficients
        cx = 1.0/(deltaX*deltaX);
        cy = 1.0/(deltaY*deltaY);
        cc = -2.0*cx - 2.0*cy - alpha;

        fX = malloc(sizeof(double) * (maxXCount));
        fY = malloc(sizeof(double) * (maxYCount));

        if (fX == NULL || fY == NULL)
        {
            fprintf(stderr, "Could not allocate memory for precalculations.");
            exit(1);
        }

        for (int x = 1; x < maxXCount-1; x++)
        {
            fX[x-1] = xStart + (x-1)*deltaX;
        }
        for (int y = 1; y < maxYCount-1; y++)
        {
            fY[y-1] = yStart + (y-1)*deltaY;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////

    double *u, *u_old;
    allocate_grid(n, m, &u, &u_old);


    double omega  = relax;
    double *src   = u_old;
    double *dst  = u;

    double error = HUGE_VAL;
    double totalError = HUGE_VAL;
    double *tmp;
    int iterationCount;
    double t1, t2;

    iterationCount = 0;
    clock_t start = clock(), diff;

    //macros are defined to translate the 2-D index (XX, YY) to the 1-dimensional array:
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]
    double updateVal;
    double f;


    char processor_name[MPI_MAX_PROCESSOR_NAME];

    ////////////////////////////////////////////////////////////////////////////////////////////////

    t1 = MPI_Wtime();


    MPI_Request RRequests[4], SRequests[4];

    MPI_Status SStatuses[4], RStatuses[4];
    /* Iterate as long as it takes to meet the convergence criterion */
    while (iterationCount < max_iteration_count && error > max_acceptable_error)
    {
        int x, y;

      /** NOTE: u(0,*), u(maxXCount-1,*), u(*,0) and u(*,maxYCount-1)
       *  are BOUNDARIES and therefore not part of the solution. Take a look at this:
       *  http://etutorials.org/Linux+systems/cluster+computing+with+linux/Part+II+Parallel+Programming/Chapter+9+Advanced+Topics+in+MPI+Programming/9.3+Revisiting+Mesh+Exchanges/
       */

        /** receive first line and column (green points) from neighbours.*/

        MPI_Irecv(&SRC(1, maxYCount-1), 1, row,    ranks.north, 0, comm_cart.id, &RRequests[0]);
        MPI_Irecv(&SRC(1, 0),           1, row,    ranks.south, 0, comm_cart.id, &RRequests[1]);
        MPI_Irecv(&SRC(maxXCount-1, 1), 1, column, ranks.east,  0, comm_cart.id, &RRequests[2]);
        MPI_Irecv(&SRC(0, 1),           1, column, ranks.west,  0, comm_cart.id, &RRequests[3]);

        /** send first line and column (green points) to neighbours.*/
        MPI_Isend(&SRC(1, 1),           1, row,    ranks.south, 0, comm_cart.id, &SRequests[0]);
        MPI_Isend(&SRC(1, maxYCount-2), 1, row,    ranks.north, 0, comm_cart.id, &SRequests[1]);
        MPI_Isend(&SRC(0, 1),           1, column, ranks.west,  0, comm_cart.id, &SRequests[2]);
        MPI_Isend(&SRC(maxXCount-2, 1), 1, column, ranks.east,  0, comm_cart.id, &SRequests[3]);

        error = 0.0;
        /** This double for loop is for white points calculations
         * Changed x and y initiate values (from 1 to 2) for white point calculations
         */
        for (y = 2; y < (maxYCount-2); y++)
        {
            for (x = 2; x < (maxXCount-2); x++)
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
        error = sqrt(error)/((maxXCount-2)*(maxYCount-2));

        MPI_Waitall(4, RRequests, RStatuses);
        /**
         * Boarder-halo points are received.
         * Calculations for green points are now made.
         * */

        //south
        y = 1;
        for (x = 1; x < maxXCount-1; x++)
        {
            updateVal = (	(SRC(x-1,y) + SRC(x+1,y))*cx +
                             (SRC(x,y-1) + SRC(x,y+1))*cy +
                             SRC(x,y)*cc
                             -(-alpha*(1.0-fX[x]*fX[x])*(1.0-fY[y]*fY[y]) - 2.0*(1.0-fX[x]*fX[x]) - 2.0*(1.0-fY[y]*fY[y])))
                        /cc;
            DST(x,y) = SRC(x,y) - omega*updateVal;
            error += updateVal*updateVal;
        }

        //north
        y = maxYCount - 2;
        for (x = 1; x < maxXCount-1; x++)
        {
            updateVal = (	(SRC(x-1,y) + SRC(x+1,y))*cx +
                             (SRC(x,y-1) + SRC(x,y+1))*cy +
                             SRC(x,y)*cc
                             -(-alpha*(1.0-fX[x]*fX[x])*(1.0-fY[y]*fY[y]) - 2.0*(1.0-fX[x]*fX[x]) - 2.0*(1.0-fY[y]*fY[y])))
                        /cc;
            DST(x,y) = SRC(x,y) - omega*updateVal;
            error += updateVal*updateVal;
        }

        //west
        x = 1;
        for (y = 1; y < maxYCount-1; y++)
        {
            updateVal = (	(SRC(x-1,y) + SRC(x+1,y))*cx +
                             (SRC(x,y-1) + SRC(x,y+1))*cy +
                             SRC(x,y)*cc
                             -(-alpha*(1.0-fX[x]*fX[x])*(1.0-fY[y]*fY[y]) - 2.0*(1.0-fX[x]*fX[x]) - 2.0*(1.0-fY[y]*fY[y])))
                        /cc;
            DST(x,y) = SRC(x,y) - omega*updateVal;
            error += updateVal*updateVal;
        }

        //east
        x = maxXCount - 2;
        for (y = 1; y < maxYCount-1; y++)
        {
            updateVal = (	(SRC(x-1,y) + SRC(x+1,y))*cx +
                             (SRC(x,y-1) + SRC(x,y+1))*cy +
                             SRC(x,y)*cc
                             -(-alpha*(1.0-fX[x]*fX[x])*(1.0-fY[y]*fY[y]) - 2.0*(1.0-fX[x]*fX[x]) - 2.0*(1.0-fY[y]*fY[y])))
                        /cc;
            DST(x,y) = SRC(x,y) - omega*updateVal;
            error += updateVal*updateVal;
        }

        //printf("\tError %g\n", error);
        iterationCount++;
        // Swap the buffers
        tmp = u_old;
        u_old = u;
        u = tmp;

        MPI_Reduce(&error, &totalError, 1, MPI_DOUBLE, MPI_SUM, 0, comm_cart.id);

        /** ALL REDUCE */
//        MPI_Allreduce(&error, &totalError, 1, MPI_DOUBLE, MPI_SUM, comm);
        MPI_Waitall(4, SRequests, SStatuses);
    }
    /** This is to avoid a WARNING about memory leaks */
    MPI_Type_free(&column);
    MPI_Type_free(&row);

    t2 = MPI_Wtime();
    printf( "Rank %d: Iterations=%3d Elapsed MPI Wall time is %f\n", comm_cart.rank, iterationCount, t2 - t1 );

    MPI_Comm_free(&comm_cart.id);
    MPI_Finalize();


    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Rank %d: Time taken %d seconds %d milliseconds\n", comm_cart.rank, msec/1000, msec%1000);

    // u_old holds the solution after the most recent buffers swap
    double absoluteError = check_solution(xLeft, yBottom,
                                          n_global+2, m_global+2,
                                          u_old,
                                          deltaX, deltaY);
    if (comm_cart.rank == 0) {
        //total error: Propably the less this value is the more accurate our solution is
        printf("Residual %g\n", totalError);
        printf("The error of the iterative solution is %g\n", absoluteError);
    }

    free(fX);
    free(fY);

    return 0;
}
