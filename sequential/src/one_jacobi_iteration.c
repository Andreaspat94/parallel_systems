#include "one_jacobi_iteration.h"
#include <stdlib.h>
#include <math.h>

#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]

void precalculate(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    double deltaX, double deltaY,
    double alpha,
    precalculations_t *pre)
{
    // Coefficients
    pre->cx = 1.0 / (deltaX*deltaX);
    pre->cy = 1.0 / (deltaY*deltaY);
    pre->cc = -2.0 * pre->cx - 2.0 * pre->cy - alpha;

    int allocCount = (maxXCount-2) * (maxYCount-2);
    pre->f = (double*)malloc(allocCount * sizeof(double));

    for (int y = 0; y < maxYCount-2; y++)
    {
        double fY = yStart + (y-1)*deltaY;
        int f_offset = y * (maxYCount - 2);

        for (int x = 0; x < maxXCount-2; x++)
        {
            double fX = xStart + (x-1)*deltaX;
            pre->f[f_offset + x] =
                - alpha*(1.0-fX*fX)*(1.0-fY*fY)
                - 2.0*(1.0-fX*fX)
                - 2.0*(1.0-fY*fY);
        }
    }
}

double one_jacobi_iteration(
    int maxXCount, int maxYCount,
    const double *src, double *dst,
    double omega,
    precalculations_t *pre)
{
    double error = 0.0;
    double updateVal;

    for (int y = 1; y < (maxYCount-1); y++)
    {
        int f_offset = (y-1) * (maxYCount - 2);

        for (int x = 1; x < (maxXCount-1); x++)
        {
            updateVal = (
                (SRC(x-1,y) + SRC(x+1,y)) * pre->cx +
                (SRC(x,y-1) + SRC(x,y+1)) * pre->cy +
                SRC(x,y) * pre->cc - pre->f[f_offset + x-1]
            ) / pre->cc;
            DST(x,y) = SRC(x,y) - omega*updateVal;
            error += updateVal*updateVal;
        }
    }

    return sqrt(error) / ((maxXCount-2)*(maxYCount-2));
}
