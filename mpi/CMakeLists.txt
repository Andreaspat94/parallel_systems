cmake_minimum_required(VERSION 3.10)
project(mpi)

set(CMAKE_C_STANDARD 11)
SET(CMAKE_C_COMPILER mpicc)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(../common/include/)

# Use all .c files in ./src/ directory and in its descendant directories.
file(GLOB_RECURSE SOURCES src/**.c ../common/src/**.c)
add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} ${MPI_C_LIBRARIES} -lm)
