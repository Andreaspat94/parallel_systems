#ifndef SEQUENTIAL_JACOBI_ITERATION_OPT3X_H
#define SEQUENTIAL_JACOBI_ITERATION_OPT3X_H

#include "./precalculations_t.h"

void jacobi_precalculate_opt3x(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    double deltaX, double deltaY,
    double alpha,
    precalculations_t *precalculations);

double jacobi_iteration_opt3x(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    const double *src, double *dst,
    double deltaX, double deltaY,
    double alpha, double omega,
    precalculations_t *precalculations);

#endif //SEQUENTIAL_JACOBI_ITERATION_OPT3X_H
