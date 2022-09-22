#include "jacobi_iteration_opt1.h"
#include <stdlib.h>
#include <math.h>

#define CX precalculations->cx
#define CY precalculations->cy
#define CC precalculations->cc
#define F  precalculations->f

void jacobi_precalculate_opt1(
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

double jacobi_iteration_opt1(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    const double *src, double *dst,
    double deltaX, double deltaY,
    double alpha, double omega,
    precalculations_t *precalculations)
{
#define FY(y) F[y]
#define FX(x) F[maxYCount-2 + x]
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]

    double error = 0.0;
    double updateVal;

    for (int y = 1; y < (maxYCount-1); y++)
    {
        double fY = FY(y-1);
        for (int x = 1; x < (maxXCount-1); x++)
        {
            double fX = FX(x-1);
            double f = -alpha*(1.0-fX*fX)*(1.0-fY*fY) - 2.0*(1.0-fX*fX) - 2.0*(1.0-fY*fY);
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
