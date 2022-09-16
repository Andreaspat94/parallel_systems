cmake_minimum_required(VERSION 3.10)
project(sequential)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_COMPILER mpicc)

# MPI is used only for measuring time.
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(../common/include/)
include_directories(./include/)

file(GLOB_RECURSE SOURCES ./src/**.c ../common/src/**.c)
message("${SOURCES}")
add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} ${MPI_C_LIBRARIES} -lm)

