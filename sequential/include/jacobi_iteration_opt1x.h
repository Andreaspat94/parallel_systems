#ifndef SEQUENTIAL_JACOBI_ITERATION_OPT1X_H
#define SEQUENTIAL_JACOBI_ITERATION_OPT1X_H

#include "./precalculations_t.h"
#include "jacobi_iteration_params_t.h"

void jacobi_precalculate_opt1x(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    double deltaX, double deltaY,
    double alpha,
    precalculations_t *precalculations);

double jacobi_iteration_opt1x(
    jacobi_iteration_params_t *params,
    precalculations_t *precalculations);

#endif //SEQUENTIAL_JACOBI_ITERATION_OPT1X_H
