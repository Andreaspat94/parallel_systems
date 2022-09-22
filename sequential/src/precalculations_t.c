#include "precalculations_t.h"
#include <stdlib.h>

void free_precalculations(precalculations_t *precalculations)
{
    free(precalculations->f);
}

