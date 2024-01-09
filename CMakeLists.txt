cmake_minimum_required(VERSION 3.24)
project(genetic_algorithm_vlsi_V2)

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -fopenmp")

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(genetic_algorithm_vlsi_V2
        main.cpp ga.cpp ga.h)
