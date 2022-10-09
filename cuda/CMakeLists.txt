cmake_minimum_required(VERSION 3.10)
project(cuda)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_COMPILER gcc)

find_package(CUDA)
include_directories(./include/)

file(GLOB_RECURSE SOURCES ./src/**.cpp ./src/**.cu ./src/**.c)
add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} m)

