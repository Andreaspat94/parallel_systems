/**
 * Precalculates all f.
 * Manual optimizations inside jacobi iteration function:
 *  As few duplicate calculations as possible (e.g. loop invariant code motion).
 *  Move the "SRC(x,y)*cc" term outside of the "updateVal" fraction; thus it becomes "SRC(x,y)".
 * 1) loop invariant code motion.
 * Memory overhead: O(n*m).
 */

#include "jacobi_iteration_opt4.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define CX precalculations->cx
#define CY precalculations->cy
#define CC precalculations->cc
#define F  precalculations->f

void jacobi_precalculate_opt4(
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
        int offset = y * maxXCount;

        for (int x = 0; x < maxXCount; x++)
        {
            double fX = xStart + x*deltaX;
            double fX2 = fX * fX;

            F[offset + x] = -alpha*(1.0-fX2)*(1.0-fY2) - 2.0*(1.0-fX2) - 2.0*(1.0-fY2);
        }
    }
}

double jacobi_iteration_opt4(
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

    for (int y = 1; y < maxYCount-1; y++)
    {
        int f_offset = (y-1) * (maxXCount-2);
        double *p_f = f + f_offset - 1; // "f_offset" and "-1" in "f[f_offset+x-1]" are loop invariants.

        int u_offset = y * maxXCount;

        for (int x = 1; x < maxXCount-1; x++)
        {
            int xx = u_offset + x;
            const double *p_src = src + xx;
            updateVal = (
                (p_src[-1]         + p_src[1])*cx +
                (p_src[-maxXCount] + p_src[maxXCount])*cy -
                p_f[x]  // previously "f[f_offset+x-1]".
            ) / cc + p_src[0];
            dst[xx] = p_src[0] - omega*updateVal;
            error += updateVal*updateVal;
        }
    }

    return sqrt(error)/((maxXCount-2)*(maxYCount-2));
}
