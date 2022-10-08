#ifndef CUDA_JACOBI_ITERATION_GPU_H
#define CUDA_JACOBI_ITERATION_GPU_H

extern "C" float jacobi_gpu(
    double *src, double *dst,
    int maxXCount, int maxYCount,
    double xStart, double yStart,
    double deltaX, double deltaY,
    double alpha, double omega,
    int max_iteration_count, double max_acceptable_error,
    int *out_iteration_count, double *out_error, float *out_elapsedTime);

#endif //CUDA_JACOBI_ITERATION_GPU_H
