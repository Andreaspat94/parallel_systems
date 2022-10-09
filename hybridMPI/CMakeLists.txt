cmake_minimum_required(VERSION 3.10)
project(openmp)

SET(CMAKE_C_STANDARD 11)
SET(CMAKE_C_COMPILER mpicc)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(../common/include/)
include_directories(./include/)

SET(CMAKE_C_COMPILER  mpicc)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Use all .c files in ./src/ directory and in its descendant directories.
file(GLOB_RECURSE SOURCES src/**.c ../common/src/**.c)
add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} PUBLIC ${MPI_C_LIBRARIES}  OpenMP::OpenMP_CXX -lm)