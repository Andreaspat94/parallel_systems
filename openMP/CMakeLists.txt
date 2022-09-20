cmake_minimum_required(VERSION 3.22)
project(parallel_systems)

set(CMAKE_C_STANDARD 11)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

SET(CMAKE_C_COMPILER  mpicc)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Use all .c files in ./src/ directory and in its descendant directories.
file(GLOB_RECURSE SOURCES src/**.c)
add_executable(parallel_systems ${SOURCES})

target_link_libraries(parallel_systems PUBLIC ${MPI_C_LIBRARIES}  OpenMP::OpenMP_CXX -lm)
target_link_libraries(parallel_systems )
