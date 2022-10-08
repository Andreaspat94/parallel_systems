#include "jacobi_gpu.h"
#include "timestamp.h"
#include <stdio.h>
#include <math.h>
#include <cuda.h>

#define THREADS_PER_BLOCK 256 // Also check with other values, e.g. 512.

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void jacobi_iteration_gpu(
    const double *d_src, double *d_dst,
    int maxXCount, int maxYCount,
    double xStart, double yStart,
    double deltaX, double deltaY,
    double alpha, double omega)
{
    int y = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int x = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int i = y * maxXCount + x;

    bool isValidIndex = y < maxYCount-1 && x < maxXCount-1;
    int ty = threadIdx.y;
    int tx = threadIdx.x;
    int ti = threadIdx.y * blockDim.x + threadIdx.x;

    __shared__ double s_src[THREADS_PER_BLOCK]; // L1 cache.

    // Each thread will do one global/L2 memory access.
    double thisCell = d_src[i];
    if (isValidIndex)
        s_src[ti] = thisCell;

    __syncthreads();

    if (!isValidIndex)
        return;

    // Corner cells will do two global/L2 reads and two L1 reads.
    // The rest border cells will do one global/L2 reads and three L1 reads.
    // All non-border cells will do four L1 reads.
    double leftCell  = tx == 0            ? d_src[i-1]         : s_src[ti-1];
    double rightCell = tx == blockDim.x-1 ? d_src[i+1]         : s_src[ti+1];
    double upperCell = ty == 0            ? d_src[i-maxXCount] : s_src[ti-blockDim.x];
    double belowCell = ty == blockDim.y-1 ? d_src[i+maxXCount] : s_src[ti+blockDim.x];

    // Coefficients
    // TODO: Is it better to always recalculate them OR precalculate them in global memory and then
    //       always accessing global memory or L2 cache?
    double cx = 1.0/(deltaX*deltaX);
    double cy = 1.0/(deltaY*deltaY);
    double cc = -2.0*cx-2.0*cy-alpha;

    // TODO: If we precalculate any of these values, then each thread will require at least one
    //       global/L2 read. Therefore, it's probably better to simply always recalculate them (?).
    double fY = yStart + y*deltaY;
    double fX = xStart + x*deltaX;
    double f = -alpha*(1.0-fX*fX)*(1.0-fY*fY) - 2.0*(1.0-fX*fX) - 2.0*(1.0-fY*fY);

    double updateVal = (
        (leftCell + rightCell) * cx +
        (upperCell + belowCell) * cy +
        thisCell * cc -f
    ) / cc;

    d_dst[i] = thisCell - omega*updateVal;

    // TODO: sum_reduce somehow all threads errors.
//    double error = 0.0;
//    error += updateVal*updateVal;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

struct pair_t {
    int first;
    int second;
};

extern "C" pair_t findMostCloseToSqrtDivisor(int number) {

    pair_t divisors = {1, number};

    for(int i = ceil(sqrt(number)); i > 0; i--)
        if(number % i == 0) {
            divisors.first = i;
            divisors.second = number/i;
            break;
        }

    return divisors;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C" float jacobi_gpu(
    double *src, double *dst,
    int maxXCount, int maxYCount,
    double xStart, double yStart,
    double deltaX, double deltaY,
    double alpha, double omega,
    int max_iteration_count, double max_acceptable_error,
    int *out_iteration_count, double *out_error, float *out_elapsedTime)
{
    double *d_src, *d_dst;
    cudaError_t err;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Allocate device memory & fill it with host memory data.

    size_t bytesCnt = maxXCount * maxYCount * sizeof(double);

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
    ///

    timestamp t_start;
    t_start = getTimestamp();

    pair_t blockSideSizes = findMostCloseToSqrtDivisor(THREADS_PER_BLOCK);
    dim3 blockSize(blockSideSizes.first, blockSideSizes.second);

    int blocksCount = ((maxXCount-2)*(maxYCount-2) + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    pair_t gridSideSizes = findMostCloseToSqrtDivisor(blocksCount);
    dim3 gridSize(gridSideSizes.first, blockSideSizes.second);

    double error = HUGE_VAL;
    int iteration_count = -1;

    while (++iteration_count < max_iteration_count && err > max_acceptable_error)
    {
        jacobi_iteration_gpu<<<gridSize, blockSize>>>(
            d_src, d_dst, maxXCount, maxYCount, xStart, yStart, deltaX, deltaY, alpha, omega);

        // TODO: calculate the total iteration's error: sqrt(error)/((maxXCount-2)*(maxYCount-2));

        // Swap buffers.
        double *tmp = d_src;
        d_src = d_dst;
        d_dst = tmp;
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
