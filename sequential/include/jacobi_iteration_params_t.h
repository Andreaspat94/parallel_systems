#ifndef SEQUENTIAL_JACOBI_ITERATION_PARAMS_T_H
#define SEQUENTIAL_JACOBI_ITERATION_PARAMS_T_H

typedef struct {
    double xStart;
    double yStart;
    int maxXCount;
    int maxYCount;
    double *src;
    double *dst;
    double deltaX;
    double deltaY;
    double alpha;
    double omega;
} jacobi_iteration_params_t;

void swap_buffers(jacobi_iteration_params_t *params);

#endif //SEQUENTIAL_JACOBI_ITERATION_PARAMS_T_H
