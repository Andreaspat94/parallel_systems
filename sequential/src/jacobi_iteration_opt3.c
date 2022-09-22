/**
 * Precalculates all f.
 * Memory overhead: O(n*m).
 */

#include "jacobi_iteration_opt3.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define CX precalculations->cx
#define CY precalculations->cy
#define CC precalculations->cc
#define F  precalculations->f

void jacobi_precalculate_opt3(
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

double jacobi_iteration_opt3(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    const double *src, double *dst,
    double deltaX, double deltaY,
    double alpha, double omega,
    precalculations_t *precalculations)
{
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]

    double error = 0.0;
    double updateVal;

    for (int y = 1; y < (maxYCount-1); y++)
    {
        int offset = (y-1) * (maxXCount-2);
        for (int x = 1; x < (maxXCount-1); x++)
        {
            double f = F[offset + x-1];
            updateVal = (
                (SRC(x-1,y) + SRC(x+1,y))*CX +
                (SRC(x,y-1) + SRC(x,y+1))*CY +
                SRC(x,y)*CC - f
            ) / CC;
            DST(x,y) = SRC(x,y) - omega*updateVal;
            error += updateVal*updateVal;
        }
    }

    return sqrt(error)/((maxXCount-2)*(maxYCount-2));
}
