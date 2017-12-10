cmake_minimum_required(VERSION 3.10)
project(BruteSHA_1)

set(CMAKE_C_STANDARD 11)

find_package(MPI)
include_directories(${MPI_C_INCLUDE_PATH})
add_definitions(-DENABLE_MPI)

add_executable(BruteSHA_1 main.c sha-1.c sha-1.h)
target_link_libraries(${PROJECT_NAME} crypto m ${MPI_C_LIBRARIES})