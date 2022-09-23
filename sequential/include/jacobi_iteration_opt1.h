#ifndef SEQUENTIAL_JACOBI_ITERATION_OPT1_H
#define SEQUENTIAL_JACOBI_ITERATION_OPT1_H

#include "./precalculations_t.h"
#include "jacobi_iteration_params_t.h"

void jacobi_precalculate_opt1(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    double deltaX, double deltaY,
    double alpha,
    precalculations_t *precalculations);

double jacobi_iteration_opt1(
    jacobi_iteration_params_t *params,
    precalculations_t *precalculations);

#endif //SEQUENTIAL_JACOBI_ITERATION_OPT1_H
