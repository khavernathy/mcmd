project(mcmd)
cmake_minimum_required(VERSION 3.5)

option(MPI "Use MPI to parallelize the calculations (requires MPI)" OFF)
option(CUDA "Use CUDA to offload force calculations to a GPU (requires CUDA)" OFF)
option(OMP "Use OpenMP to parallelize the calculations (requires OpenMP)" OFF)
option(WINDOWS "Install MCMD on Windows" OFF)

add_definitions( -D`echo VERSION=\\`git rev-list HEAD|wc -l\\``)

configure_file (
"${PROJECT_SOURCE_DIR}/src/cmake_config.h.in"
"${PROJECT_BINARY_DIR}/src/cmake_config.h"
)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "-lm -I. -Ofast -foptimize-sibling-calls -finline-limit=10000 -fexpensive-optimizations -flto -frename-registers")

set(CMAKE_CXX_STANDARD 11)

set(INCLUDE src/)

set(SRC 
src/main.cpp
)


if(MPI)
    message("-- MPI Enabled")
    find_package(MPI REQUIRED)
    if(NOT MPI_C_FOUND)
        message(FATAL_ERROR "-- MPI not found! Exiting ...")
    endif()
    set(INCLUDE ${INCLUDE} ${MPI_C_INCLUDE_PATH})
    set(LIB ${LIB} ${MPI_C_LIBRARIES})
else()
    message("-- MPI Disabled")
endif()


if(OMP)
    message("-- OpenMP Enabled")
    find_package(OPENMP REQUIRED)
    if(NOT OPENMP_C_FOUND)
        message(FATAL_ERROR "-- OpenMP not found! Exiting ...")
    endif()
else()
    message("-- OpenMP Disabled")
endif()

if(CUDA)
    message("-- CUDA Enabled")
    find_package(CUDA REQUIRED)
    set(SRC ${SRC} src/polarization_gpu/polar_cuda.cu)
else()
    message("-- CUDA Disabled")
endif()




# make executable

if(CUDA)
    cuda_add_executable(${PROJECT_NAME} ${SRC})
else()
    add_executable(${PROJECT_NAME} ${SRC})
endif()

