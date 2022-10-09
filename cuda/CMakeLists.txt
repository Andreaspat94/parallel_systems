cmake_minimum_required(VERSION 3.10)
project(cuda)

set(CMAKE_CXX_STANDARD 11)

include_directories(../common/include/)
include_directories(./include/)

file(GLOB_RECURSE SOURCES ./src/**.cpp ./src/**.cu ../common/src/**.c)
add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} m)

