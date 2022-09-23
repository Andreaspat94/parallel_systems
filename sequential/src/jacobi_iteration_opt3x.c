/**
 * Precalculates all "f".
 * Memory overhead: O(n*m).
 *
 * Manual optimizations inside jacobi iteration function:
 * - As few duplicate calculations as possible (e.g. loop invariant code motion).
 * - Move the "SRC(x,y)*cc" term outside of the "updateVal" fraction; thus it becomes "SRC(x,y)".
 */

#include "jacobi_iteration_opt3x.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define CX precalculations->cx
#define CY precalculations->cy
#define CC precalculations->cc
#define F  precalculations->f

void jacobi_precalculate_opt3x(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    double deltaX, double deltaY,
    double alpha,
    precalculations_t *precalculations)
{
    // Coefficients
    CX = 1.0/(deltaX*deltaX);
    CY = 1.0/(deltaY*deltaY);
    CC = -2.0*CX-2.0*CY-alpha;

    F = malloc(sizeof(double) * (maxYCount * maxXCount));
    if (F == NULL)
    {
        fprintf(stderr, "Could not allocate memory for precalculations.");
        exit(1);
    }

    for (int y = 0; y < maxYCount; y++)
    {
        double fY = yStart + y*deltaY;
        double fY2 = fY * fY;
        int x_offset = y * maxXCount;

        for (int x = 0; x < maxXCount; x++)
        {
            double fX = xStart + x*deltaX;
            double fX2 = fX * fX;

            F[x_offset + x] = -alpha*(1.0-fX2)*(1.0-fY2) - 2.0*(1.0-fX2) - 2.0*(1.0-fY2);
        }
    }
}

double jacobi_iteration_opt3x(
    jacobi_iteration_params_t *params,
    precalculations_t *precalculations)
{
    int maxXCount = params->maxXCount;
    int maxYCount = params->maxYCount;
    double *src   = params->src;
    double *dst   = params->dst;
    double omega  = params->omega;

    double error = 0.0;
    double updateVal;
    double cx = CX;
    double cy = CY;
    double cc = CC;
    double *f = F;

    for (int y = 1; y < maxYCount-1; y++)
    {
        int x_offset = y * maxXCount;
        double *_f = f + (y-1)*(maxXCount-2) - 1; // The "-1" simplifies "_f[x-1]" to "_f[x]".

        for (int x = 1; x < maxXCount-1; x++)
        {
            int _x = x_offset + x;
            const double *_src = src + _x;

            updateVal = (
                (_src[-1]         + _src[1])*cx +
                (_src[-maxXCount] + _src[maxXCount])*cy -
                _f[x]
            ) / cc + _src[0];
            dst[_x] = _src[0] - omega*updateVal;
            error += updateVal*updateVal;
        }
    }

    return sqrt(error)/((maxXCount-2)*(maxYCount-2));
}
