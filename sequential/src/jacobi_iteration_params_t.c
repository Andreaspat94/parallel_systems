#include "jacobi_iteration_params_t.h"

void swap_buffers(jacobi_iteration_params_t *params)
{
    double *tmp = params->src;
    params->src = params->dst;
    params->dst = tmp;
}

