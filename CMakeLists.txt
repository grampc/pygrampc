cmake_minimum_required(VERSION 3.15...3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX C)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Eigen 3.4 is required
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

include_directories(
    include(
        src/
        grampc/include
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
    grampc/src/trapezodial.c
)

set(PYGRAMPC_SOURCES
    src/problem_description.cpp
    src/grampc_interface.cpp
)

pybind11_add_module(_core MODULE ${GRAMPC_SOURCES} ${PYGRAMPC_SOURCES} src/binding.cpp)
target_link_libraries(_core PRIVATE Eigen3::Eigen)

install(TARGETS _core DESTINATION pygrampc)
