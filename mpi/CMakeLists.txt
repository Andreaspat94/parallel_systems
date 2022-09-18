cmake_minimum_required(VERSION 3.22)
project(parallel_systems)

set(CMAKE_C_STANDARD 11)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

SET(CMAKE_C_COMPILER  mpicc)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# Use all .c files in ./src/ directory and in its descendant directories.
file(GLOB_RECURSE SOURCES src/**.c)
add_executable(parallel_systems ${SOURCES})

target_link_libraries(parallel_systems ${MPI_C_LIBRARIES} -lm)
