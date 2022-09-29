#ifndef SEQUENTIAL_PRINTS_H
#define SEQUENTIAL_PRINTS_H

#include "input.h"
#include "timing.h"

void print_input(const input_t *input);
void print_iterations(int iteration_count);
void print_times(const times_t *times);
void print_residual(double error);
void print_iterative_solution_error(double absolute_error);
void print_output(int iteration_count, const times_t *times, double error, double absolute_error);

#endif //SEQUENTIAL_PRINTS_H
