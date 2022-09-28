#include "../include/common/allocate_grid.h"
#include <stdlib.h>
#include <stdio.h>

void allocate_grid(int maxXCount, int maxYCount, double **u, double **u_old)
{
    int allocCount = maxXCount * maxYCount;

    // Those two calls also zero the boundary elements
    *u     = (double*)calloc(allocCount, sizeof(double)); //reverse order
    *u_old = (double*)calloc(allocCount, sizeof(double));

    // printf("allocCount=%d u=%p u_old=%p\n", allocCount, u, u_old);

    if (u == NULL || u_old == NULL)
    {
        printf("Not enough memory for two %ix%i matrices\n", maxXCount, maxYCount);
        exit(1);
    }
}
