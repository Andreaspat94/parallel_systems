#ifndef SEQUENTIAL_ONE_JACOBI_ITERATION_H
#define SEQUENTIAL_ONE_JACOBI_ITERATION_H

typedef struct {
    double cx;
    double cy;
    double cc;
    double *f;
} precalculations_t;

void precalculate(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    double deltaX, double deltaY,
    double alpha,
    precalculations_t *pre);

/**
 * Performs one iteration of the Jacobi method and computes the residual value.
 *
 * NOTE: u(0,*), u(maxXCount-1,*), u(*,0) and u(*,maxYCount-1) are BOUNDARIES and therefore not part
 * of the solution.
 */
double one_jacobi_iteration(
    int maxXCount, int maxYCount,
    const double *src, double *dst,
    double omega,
    precalculations_t *pre);

#endif //SEQUENTIAL_ONE_JACOBI_ITERATION_H
