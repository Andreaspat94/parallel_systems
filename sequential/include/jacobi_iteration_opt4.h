#ifndef SEQUENTIAL_JACOBI_ITERATION_OPT4_H
#define SEQUENTIAL_JACOBI_ITERATION_OPT4_H

#include "./precalculations_t.h"

void jacobi_precalculate_opt4(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    double deltaX, double deltaY,
    double alpha,
    precalculations_t *precalculations);

double jacobi_iteration_opt4(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    const double *src, double *dst,
    double deltaX, double deltaY,
    double alpha, double omega,
    precalculations_t *precalculations);

#endif //SEQUENTIAL_JACOBI_ITERATION_OPT4_H
