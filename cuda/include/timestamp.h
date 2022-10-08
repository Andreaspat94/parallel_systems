#ifndef CUDA_TIMESTAMP_H_
#define CUDA_TIMESTAMP_H_

#ifdef __linux__

#include <cstdlib>
#include <sys/time.h>

typedef timeval timestamp;

inline timestamp getTimestamp() {
    timeval t{};
    gettimeofday(&t, nullptr);
    return t;
}

inline float getElapsedtime(timestamp t) {
    timeval tn{};
    gettimeofday(&tn, nullptr);
    return (float)(tn.tv_sec - t.tv_sec) * 1000.0f + (float)(tn.tv_usec - t.tv_usec) / 1000.0f;
}

#else

#include <ctime>

typedef clock_t timestamp;

inline timestamp getTimestamp(void) {
    return clock();
}

inline float getElapsedtime(timestamp t) {
    return ((float) clock() - t) / CLOCKS_PER_SEC * 1000.0f;
}

#endif

#endif //CUDA_TIMESTAMP_H_
