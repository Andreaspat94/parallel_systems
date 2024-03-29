/**
 * Precalculates all "fX^2" and "fY^2".
 * Memory overhead: O(n+m).
 */

#include "jacobi_iteration_opt2.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define CX precalculations->cx
#define CY precalculations->cy
#define CC precalculations->cc
#define F  precalculations->f

void jacobi_precalculate_opt2(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    double deltaX, double deltaY,
    double alpha,
    precalculations_t *precalculations)
{
#define FY2(y) F[y]
#define FX2(x) F[maxYCount + x]

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
        double fY = yStart + y*deltaY;
        FY2(y) = fY * fY;
    }

    for (int x = 0; x < maxXCount; x++)
    {
        double fX = xStart + x*deltaX;
        FX2(x) = fX * fX;
    }

#undef FY2
#undef FX2
}

double jacobi_iteration_opt2(
    jacobi_iteration_params_t *params,
    precalculations_t *precalculations)
{
#define FY2(y) F[y]
#define FX2(x) F[maxYCount-2 + x]
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]

    int maxXCount = params->maxXCount;
    int maxYCount = params->maxYCount;
    double *src   = params->src;
    double *dst   = params->dst;
    double alpha  = params->alpha;
    double omega  = params->omega;

    double error = 0.0;
    double updateVal;

    for (int y = 1; y < (maxYCount-1); y++)
    {
        double fY2 = FY2(y-1);
        for (int x = 1; x < (maxXCount-1); x++)
        {
            double fX2 = FX2(x-1);
            double f = -alpha*(1.0-fX2)*(1.0-fY2) - 2.0*(1.0-fX2) - 2.0*(1.0-fY2);
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
