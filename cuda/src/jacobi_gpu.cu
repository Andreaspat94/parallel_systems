#include "jacobi_gpu.h"
#include "timestamp.h"
#include <math.h>
#include <cuda.h>

extern "C" double jacobi_iteration_gpu(
    const double *src, double *dst,
    int maxXCount, int maxYCount,
    double xStart, double yStart,
    double deltaX, double deltaY,
    double alpha, double omega)
{
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]

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

extern "C" float jacobi_gpu(
    double *src, double *dst,
    int maxXCount, int maxYCount,
    double xStart, double yStart,
    double deltaX, double deltaY,
    double alpha, double omega,
    int max_iteration_count, double max_acceptable_error,
    int *out_iteration_count, double *out_error, float *out_elapsedTime)
{
    float *d_src, *d_dst;
    cudaError_t err;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Allocate device memory & fill it with host memory data.

    size_t bytesCnt = maxXCount * maxYCount * sizeof(float);

    // Allocate device memory.
    err = cudaMalloc((void **) &d_src, bytesCnt);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s\n", err);
        return err;
    }
    err = cudaMalloc((void **) &d_dst, bytesCnt);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s\n", err);
        return err;
    }

    // Copy data to device memory.
    err = cudaMemcpy(d_src, src, bytesCnt, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s\n", err);
        return err;
    }
    err = cudaMemcpy(d_dst, dst, bytesCnt, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s\n", err);
        return err;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Run the jacobi iterations.

    double error = HUGE_VAL;
    int iteration_count = -1;

    timestamp t_start;
    t_start = getTimestamp();

    while (++iteration_count < max_iteration_count && err > max_acceptable_error)
    {
        error = jacobi_iteration_gpu(
            src, dst, maxXCount, maxYCount, xStart, yStart, deltaX, deltaY, alpha, omega);

        // Swap buffers.
        double *temp = src; src = dst; dst = temp;
    }

    float elapsedTime = getElapsedtime(t_start);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Copy result back to host memory & deallocate device memory.

    // Copy results back to host memory
    err = cudaMemcpy(dst, d_dst, bytesCnt, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s\n", err);
        return err;
    }

    err = cudaFree(d_dst);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s\n", err);
        return err;
    }
    err = cudaFree(d_src);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s\n", err);
        return err;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Set function output values and return ok.

    *out_elapsedTime = elapsedTime;
    *out_iteration_count = iteration_count;
    *out_error = error;

    return 0.f;
}
