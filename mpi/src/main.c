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
#include <string.h>
#include "common/read_input.h"
#include "common/allocate_grid.h"

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
} neighbour_ranks_t;

bool is_perfect_square(int number)
{
    double root = sqrt(number);
    return root == floor(root);
}

int main(int argc, char **argv)
{
    bool print_all_process_times = false;

    for (int i = 1; i < argc; i++)
    {
        if (!strcmp("--all-process-times", argv[i]))
        {
            print_all_process_times = true;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Initialize MPI and collect MPI_COMM_WORLD-related info.

    MPI_Init(NULL, NULL);

    comm_t comm_world = { MPI_COMM_WORLD };
    MPI_Comm_size(comm_world.id, &comm_world.size);
    MPI_Comm_rank(comm_world.id, &comm_world.rank);

    // Require that the total number of processes ∈ {1,4,9,16,25,36,49,64,80}.
    if (comm_world.size != 80 && !(is_perfect_square(comm_world.size) && comm_world.size < 80))
    {
        MPI_Abort(comm_world.id, 1);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Read input values from stdin (only main process) and broadcast them to all processes.

    int n_global, m_global;
    double alpha, omega;
    double max_acceptable_error;
    int max_iteration_count;

    bcast_input_t input;

    if (comm_world.rank == 0)
    {
        read_input(&input.n, &input.m, &input.alpha, &input.relax,
                   &input.max_acceptable_error, &input.max_iteration_count, true);
    }

    MPI_Bcast(&input, sizeof(bcast_input_t), MPI_BYTE, 0, comm_world.id);

    n_global = input.n;
    m_global = input.m;
    alpha = input.alpha;
    omega = input.relax;
    max_acceptable_error = input.max_acceptable_error;
    max_iteration_count = input.max_iteration_count;

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

    neighbour_ranks_t ranks; // Stencil/neighbour process ranks inside the topology.

    // direction = 0 or 1 corresponding to the two dimensions x,y.
    MPI_Cart_shift(comm_cart.id, 0, 1, &ranks.west,  &ranks.east);
    MPI_Cart_shift(comm_cart.id, 1, 1, &ranks.north, &ranks.south);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Create the grids "u" and "u_old" and relevant MPI datatypes for managing their rows and
    /// columns.

    // (n x m) is the per-process size of its data sub-grid.
    int n = n_global / dims[0];
    int m = m_global / dims[1];

    // (maxXCount x maxYCount) is the per-process size of its data sub-grid including the halos.
    int maxXCount = n + 2;
    int maxYCount = m + 2;

    double *u, *u_old;
    allocate_grid(maxXCount, maxYCount, &u, &u_old);

    MPI_Datatype row;
    MPI_Type_contiguous(n, MPI_DOUBLE, &row);
    MPI_Type_commit(&row);

    MPI_Datatype column;
    MPI_Type_vector(m, 1, maxXCount, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Solve in [-1, 1] x [-1, 1].

    double xLeft   = -1.0, xRight = 1.0;
    double yBottom = -1.0,    yUp = 1.0;

    double xStart = xLeft;
    double yStart = yBottom;

    double deltaX = (xRight-xLeft)/(n_global-1);
    double deltaY = (yUp-yBottom)/(m_global-1);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Make any jacobi-iteration precalculations.

    // Coefficients
    double cx = 1.0/(deltaX*deltaX);
    double cy = 1.0/(deltaY*deltaY);
    double cc = -2.0*cx - 2.0*cy - alpha;

    double *fX = malloc(sizeof(double) * maxXCount);
    double *fY = malloc(sizeof(double) * maxYCount);

    if (fX == NULL || fY == NULL)
    {
        fprintf(stderr, "Could not allocate memory for precalculations.");
        MPI_Abort(comm_cart.id, 1);
    }

    for (int x = 0; x < n; x++)
    {
        fX[x+1] = xStart + (coords[0]*n + x)*deltaX;
    }
    for (int y = 0; y < m; y++)
    {
        fY[y+1] = yStart + (coords[1]*m + y)*deltaY;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Start executing the jacobi iterations.

    double *src = u_old;
    double *dst = u;
    double *tmp;
    double update_val;
    double error_global = HUGE_VAL;
    int iteration_count = 0;

    MPI_Barrier(comm_cart.id);

    double t1 = MPI_Wtime();
    clock_t clock1 = clock();

    MPI_Request recv_requests[4], send_requests[4];
    MPI_Status send_statuses[4];

    // Calculate the x and y ranges that the double for loop will operate upon.
    //  NOTE: A negative rank means that the process has no neighbour at that specific side.
    int wp_y_begin = ranks.north < 0 ? 1 : 2;
    int wp_y_end   = ranks.south < 0 ? maxYCount-1 : maxYCount-2;
    int wp_x_begin = ranks.west < 0  ? 1 : 2;
    int wp_x_end   = ranks.east < 0  ? maxXCount-1 : maxXCount-2;

    while (iteration_count < max_iteration_count && error_global > max_acceptable_error)
    {
        // These two macros translate the 2-D index (XX, YY) to the 1-dimensional array:
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]

        // NOTE: u(0,*), u(maxXCount-1,*), u(*,0) and u(*,maxYCount-1) are BOUNDARIES and therefore
        // not part of the solution. Take a look at this:
        // http://etutorials.org/Linux+systems/cluster+computing+with+linux/Part+II+Parallel+Programming/Chapter+9+Advanced+Topics+in+MPI+Programming/9.3+Revisiting+Mesh+Exchanges/

        // Receive adjacent lines and columns from neighbours to fill my halo points.
        MPI_Irecv(&SRC(1, 0),           1, row,    ranks.north, 0, comm_cart.id, &recv_requests[0]);
        MPI_Irecv(&SRC(1, maxYCount-1), 1, row,    ranks.south, 0, comm_cart.id, &recv_requests[1]);
        MPI_Irecv(&SRC(0, 1),           1, column, ranks.west,  0, comm_cart.id, &recv_requests[2]);
        MPI_Irecv(&SRC(maxXCount-1, 1), 1, column, ranks.east,  0, comm_cart.id, &recv_requests[3]);

        // Send my border lines and columns to neighbours.
        MPI_Isend(&SRC(1, 1),           1, row,    ranks.north, 0, comm_cart.id, &send_requests[0]);
        MPI_Isend(&SRC(1, maxYCount-2), 1, row,    ranks.south, 0, comm_cart.id, &send_requests[1]);
        MPI_Isend(&SRC(1, 1),           1, column, ranks.west,  0, comm_cart.id, &send_requests[2]);
        MPI_Isend(&SRC(maxXCount-2, 1), 1, column, ranks.east,  0, comm_cart.id, &send_requests[3]);

#define UPDATE_VAL(XX,YY) ((\
    (SRC((XX)-1,(YY)) + SRC((XX)+1,(YY)))*cx +\
    (SRC((XX),(YY)-1) + SRC((XX),(YY)+1))*cy +\
    SRC((XX),(YY))*cc\
    -(-alpha*(1.0-fX[XX]*fX[XX])*(1.0-fY[YY]*fY[YY]) - 2.0*(1.0-fX[XX]*fX[XX]) - 2.0*(1.0-fY[YY]*fY[YY]))\
) / cc)

        double error = 0.0;

        // Calculate all white points.
        // Also calculate the green points on the sides that there are no neighbours.
        for (int y = wp_y_begin; y < wp_y_end; y++)
        {
            for (int x = wp_x_begin; x < wp_x_end; x++)
            {
                update_val = UPDATE_VAL(x,y);
                DST(x,y) = SRC(x,y) - omega*update_val;
                error += update_val*update_val;
            }
        }

        // Calculate the green points on the sides that there are neighbours.
        // The for loop's logic is the following:
        //  1. Wait for any neighbour process to finish sending the data we need.
        //  2. Calculate the appropriate green points based on the data we received.
        //  3. If there are remaining neighbours that have not yet sent their data, then go to
        //     step 1; otherwise this loop is completed.
        for (int i = 1; i <= 4; i++)
        {
            int index;
            MPI_Status status;

            MPI_Waitany(4, recv_requests, &index, &status);

            if (status.MPI_SOURCE < 0)
                continue;

            bool flag = false;

            if ((flag = status.MPI_SOURCE == ranks.north) || status.MPI_SOURCE == ranks.south)
            {
                int y = flag ? 1 : maxYCount-2; // Top or bottom row.
                for (int x = 1; x < maxXCount-1; x++)
                {
                    update_val = UPDATE_VAL(x,y);
                    DST(x,y) = SRC(x,y) - omega*update_val;
                    error += update_val*update_val;
                }
            }
            else if ((flag = status.MPI_SOURCE == ranks.west) || status.MPI_SOURCE == ranks.east)
            {
                int x = flag ? 1 : maxXCount-2; // Left or right column.
                for (int y = 1; y < maxYCount-1; y++)
                {
                    update_val = UPDATE_VAL(x,y);
                    DST(x,y) = SRC(x,y) - omega*update_val;
                    error += update_val*update_val;
                }
            }
        }

        // Calculate the iteration's total error.
        double error_iteration_sum;
        MPI_Allreduce(&error, &error_iteration_sum, 1, MPI_DOUBLE, MPI_SUM, comm_cart.id);
        error_global = sqrt(error_iteration_sum)/(n_global*m_global);

        //if (comm_cart.rank == 0)
        //    printf("===============> %g\n", error_global);

        iteration_count++;

        // Swap the buffers
        tmp = src; src = dst; dst = tmp;

        MPI_Waitall(4, send_requests, send_statuses); // TODO: do we need this? Maybe not...
    }

    double t2 = MPI_Wtime();
    clock_t clock2 = clock();

    double wtime_dt = t2 - t1;
    clock_t clock_dt = (clock2 - clock1) * 1000 / CLOCKS_PER_SEC;

    free(fX);
    free(fY);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Print results.

    if (comm_cart.rank == 0)
    {
        printf("-> Residual %g\n", error_global);
    }

    if (print_all_process_times)
    {
        printf("-> [rank=%2d] Iterations: %2d, MPI_Wtime: %f secs, clock: %ld.%03ld secs\n",
               comm_cart.rank, iteration_count, wtime_dt, clock_dt/1000, clock_dt%1000);
    }
    else
    {
        double wtime_dt_max;
        MPI_Reduce(&wtime_dt, &wtime_dt_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm_cart.id);
        if (comm_cart.rank == 0)
        {
            printf("-> Iterations: %2d, MPI_Wtime: %f secs, clock: %ld.%03ld secs\n",
                   iteration_count, wtime_dt, clock_dt/1000, clock_dt%1000);
        }
    }

    MPI_Type_free(&row);
    MPI_Type_free(&column);
    MPI_Comm_free(&comm_cart.id);

    // TODO: also parallelize check_solution?
//    // u_old holds the solution after the most recent buffers swap
//    double absoluteError = check_solution(xLeft, yBottom,
//                                          n_global+2, m_global+2,
//                                          u_old,
//                                          deltaX, deltaY);
//
//    if (comm_world.rank == 0) {
//        //total error: Propably the less this value is the more accurate our solution is
//        printf("Residual %g\n", totalError);
//        printf("The error of the iterative solution is %g\n", absoluteError);
//    }

    free(u);
    free(u_old);

    MPI_Finalize();

    return 0;
}
