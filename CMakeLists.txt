cmake_minimum_required(VERSION 3.24)
project(genetic_algorithm_vlsi_V2)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(genetic_algorithm_vlsi_V2
        main.cpp fitness.cpp fitness.h population.cpp population.h ga.cpp ga.h)
