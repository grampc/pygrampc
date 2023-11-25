cmake_minimum_required(VERSION 3.15...3.26)
project(my_problem LANGUAGES CXX C)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

include_directories(
    include(
        ../../eigen/
        ../../src/
        ../../src/GRAMPC/include
    )
)

pybind11_add_module(my_problem MODULE MyProblem.cpp)
