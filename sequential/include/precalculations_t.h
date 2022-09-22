#ifndef SEQUENTIAL_PRECALCULATIONS_T_H
#define SEQUENTIAL_PRECALCULATIONS_T_H

typedef struct {
    double cx;
    double cy;
    double cc;
    double *f;
} precalculations_t;

void free_precalculations(precalculations_t *precalculations);

#endif //SEQUENTIAL_PRECALCULATIONS_T_H
