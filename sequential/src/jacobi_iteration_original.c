#include "jacobi_iteration_original.h"
#include <math.h>
#include <stdlib.h>

void jacobi_precalculate_original(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    double deltaX, double deltaY,
    double alpha,
    precalculations_t *precalculations)
{
    precalculations->f = NULL;
}

double jacobi_iteration_original(
    jacobi_iteration_params_t *params,
    precalculations_t *precalculations)
{
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]

    double xStart = params->xStart;
    double yStart = params->yStart;
    int maxXCount = params->maxXCount;
    int maxYCount = params->maxYCount;
    double *src   = params->src;
    double *dst   = params->dst;
    double deltaX = params->deltaX;
    double deltaY = params->deltaY;
    double alpha  = params->alpha;
    double omega  = params->omega;

    // Coefficients
    double cx = 1.0/(deltaX*deltaX);
    double cy = 1.0/(deltaY*deltaY);
    double cc = -2.0*cx-2.0*cy-alpha;

    double error = 0.0;
    double updateVal;

    for (int y = 1; y < (maxYCount-1); y++)
    {
        double fY = yStart + (y-1)*deltaY;
        for (int x = 1; x < (maxXCount-1); x++)
        {
            double fX = xStart + (x-1)*deltaX;
            double f = -alpha*(1.0-fX*fX)*(1.0-fY*fY) - 2.0*(1.0-fX*fX) - 2.0*(1.0-fY*fY);
            updateVal = (
                (SRC(x-1,y) + SRC(x+1,y))*cx +
                (SRC(x,y-1) + SRC(x,y+1))*cy +
                SRC(x,y)*cc - f
            ) / cc;
            DST(x,y) = SRC(x,y) - omega*updateVal;
            error += updateVal*updateVal;
        }
    }

    return sqrt(error)/((maxXCount-2)*(maxYCount-2));
}
