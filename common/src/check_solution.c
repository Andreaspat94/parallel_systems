#include "../include/common/check_solution.h"

#define U(XX,YY) u[(YY)*maxXCount+(XX)]

double check_solution(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    const double *u,
    double deltaX, double deltaY)
{
    int x, y;
    double fX, fY;
    double localError, error = 0.0;

    for (y = 1; y < (maxYCount-1); y++)
    {
        fY = yStart + (y-1)*deltaY;

        for (x = 1; x < (maxXCount-1); x++)
        {
            fX = xStart + (x-1)*deltaX;
            localError = U(x,y) - (1.0-fX*fX)*(1.0-fY*fY);
            error += localError*localError;
        }
    }
    return error;
}
