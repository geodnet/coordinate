cmake_minimum_required(VERSION 3.16)
project(coordinate)

set(CMAKE_CXX_STANDARD 17)

add_executable(coordinate
coordinate.cpp
coord.h
coord.c
coord_sol.h
coord_sol.cpp
    )

target_link_libraries(coordinate m)