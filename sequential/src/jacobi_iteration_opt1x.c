/**
 * Precalculates all "fX" and "fY".
 * Memory overhead: O(n+m).
 *
 * Manual optimizations inside jacobi iteration function:
 * - As few duplicate calculations as possible (e.g. loop invariant code motion).
 * - Move the "SRC(x,y)*cc" term outside of the "updateVal" fraction; thus it becomes "SRC(x,y)".
 * - Simplify calculation of "f" (find an equal but more performant formula).
 */

#include "jacobi_iteration_opt1x.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define CX precalculations->cx
#define CY precalculations->cy
#define CC precalculations->cc
#define F  precalculations->f

void jacobi_precalculate_opt1x(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    double deltaX, double deltaY,
    double alpha,
    precalculations_t *precalculations)
{
#define FY(y) F[y]
#define FX(x) F[maxYCount + x]

    // Coefficients
    CX = 1.0/(deltaX*deltaX);
    CY = 1.0/(deltaY*deltaY);
    CC = -2.0*CX-2.0*CY-alpha;

    F = malloc(sizeof(double) * (maxYCount + maxXCount));
    if (F == NULL)
    {
        fprintf(stderr, "Could not allocate memory for precalculations.");
        exit(1);
    }

    for (int y = 0; y < maxYCount; y++)
    {
        FY(y) = yStart + y*deltaY;
    }

    for (int x = 0; x < maxXCount; x++)
    {
        FX(x) = xStart + x*deltaX;
    }

#undef FY
#undef FX
}

double jacobi_iteration_opt1x(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    const double *src, double *dst,
    double deltaX, double deltaY,
    double alpha, double omega,
    precalculations_t *precalculations)
{
    double error = 0.0;
    double updateVal;
    double cx = CX;
    double cy = CY;
    double cc = CC;
    double *f = F;
    double alpha_plus_2 = alpha + 2;
    double alpha_plus_4 = alpha + 4;

    for (int y = 1; y < maxYCount-1; y++)
    {
        int x_offset = y * maxXCount;
        double fY = f[y-1];
        double *_f = f + maxYCount - 2;

        for (int x = 1; x < maxXCount-1; x++)
        {
            int _x = x_offset + x;
            const double *_src = src + _x;
            double fX = _f[x-1];
            double fv = alpha_plus_2*(fX*fX + fY*fY) - alpha*fX*fX*fY*fY - alpha_plus_4;

            updateVal = (
                (_src[-1]         + _src[1])*cx +
                (_src[-maxXCount] + _src[maxXCount])*cy -
                fv
            ) / cc + _src[0];
            dst[_x] = _src[0] - omega*updateVal;
            error += updateVal*updateVal;
        }
    }

    return sqrt(error)/((maxXCount-2)*(maxYCount-2));
}
