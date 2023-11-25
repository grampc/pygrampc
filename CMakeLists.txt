cmake_minimum_required(VERSION 3.15...3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX C)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

include_directories(
    include(
        src/
        grampc/include
        eigen/
    )
)

set(GRAMPC_SOURCES
    grampc/src/euler1.c
    grampc/src/eulermod2.c
    grampc/src/grampc_alloc.c
    grampc/src/grampc_fixedsize.c
    grampc/src/grampc_init.c
    grampc/src/grampc_mess.c
    grampc/src/grampc_run.c
    grampc/src/grampc_setopt.c
    grampc/src/grampc_setparam.c
    grampc/src/grampc_util.c
    grampc/src/heun2.c
    grampc/src/rodas.c
    grampc/src/ruku45.c
    grampc/src/simpson.c
    grampc/src/trapezodial.c)

pybind11_add_module(_core MODULE ${GRAMPC_SOURCES} src/problem_base.cpp src/grampc_binding.cpp)

install(TARGETS _core DESTINATION pygrampc)
