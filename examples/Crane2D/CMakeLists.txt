cmake_minimum_required(VERSION 3.15...3.26)
project(crane_problem LANGUAGES CXX C)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

include_directories(
    include(
        ../../eigen/
        ../../src/
        ../../grampc/include
    )
)

pybind11_add_module(crane_problem MODULE Crane2D.cpp)
