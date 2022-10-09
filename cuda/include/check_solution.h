#ifndef SEQUENTIAL_CHECK_SOLUTION_H
#define SEQUENTIAL_CHECK_SOLUTION_H

/**
 * Checks the error between numerical and exact solutions.
 */
double check_solution(
    double xStart, double yStart,
    int maxXCount, int maxYCount,
    const double *u,
    double deltaX, double deltaY);

#endif //SEQUENTIAL_CHECK_SOLUTION_H
