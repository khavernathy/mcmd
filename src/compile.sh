#!/bin/bash

#########################
# HOW TO USE THIS SCRIPT
# bash compile.sh [type] [HPC cluster if any]
#
# e.g.
#
# bash compile.sh               Defaults to cpu compilation
# bash compile.sh cpu           (same result; for a single computer on Mac or Linux)
# bash compile.sh cpu linux     (optimized for linux (ONLY))
# bash compile.sh cpu windows   (for use on Windows ONLY -- you must have gcc installed, e.g. through Cygwin)
# bash compile.sh cpu debug     (same but with errors/warnings)
# bash compile.sh cpu circe     (for CIRCE (HPC at USF))
# bash compile.sh cpu bridges   (for bridges (HPC at UPitt)

# bash compile.sh icpu bridges  (for bridges Intel compilation (seems slower than gcc))

# bash compile.sh gpu           (for a single computer with CUDA compatible GPU installed).
# bash compile.sh gpu debug     (for GPU compilation with errors)
# bash compile.sh gpu circe     (CUDA GPU on circe HPC)

# bash compile.sh omp           (including OpenMP support)
# bash compile.sh omp linux     (same, optimized for Linux)
##########################
###echo "This should take like 10 sec."


DEFAULT="cpu"

if [ $# -eq 0 ]; then # IF NO ARGUMENT GIVEN
    option=$DEFAULT   # MANUAL OPTION (only reads if no argument given)
elif [ $1 == "cpu" ]; then
    option="cpu"
elif [ $1 == "gpu" ]; then
    option="gpu";
elif [ $1 == "mpi" ]; then
    option="mpi"
elif [ $1 == "omp" ]; then
    option="omp"
else
echo "Invalid options detected. Not compiling. Read the header comments for compilation examples."
fi

if [[ "$option" == "cpu" ]]; then
    # THIS IS FOR SERIAL COMPILATION (1 CPU ONLY, NO GPU)
    echo "Doing serial GCC (1 processor) compilation for CPU"
    if [[ "$2" == "circe" ]]; then
        echo "... for CIRCE cluster environment.";
        module purge
        module load compilers/gcc/6.2.0
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -foptimize-sibling-calls -finline-limit=10000 -fexpensive-optimizations -flto -march=native -frename-registers
    elif [[ "$2" == "bridges" ]]; then
        echo "... for Bridges cluster environment.";
        module purge
        module load gcc/6.3.0
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -foptimize-sibling-calls -finline-limit=10000 -fexpensive-optimizations -flto -march=native -frename-registers
    elif [[ "$2" == "debug" ]]; then
        echo "... in debug mode (with errors/warnings";
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -Werror -Wall;
    elif [[ "$2" == "linux" ]]; then
        echo "... optimized for linux."
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -foptimize-sibling-calls -finline-limit=10000 -fexpensive-optimizations -flto -march=native -frename-registers 
    elif [[ "$2" == "O3" ]]; then
        echo "... optimized at -O3 (not fully)."
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -O3
    elif [[ "$2" == "windows" ]]; then
        echo "... for Windows OS.";
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -D WINDOWS
    else
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast;
    fi
######################################################################
elif [[ "$option" == "gpu" ]]; then
    # GPU compilation enabled
    echo "Doing serial GCC (1 processor) compilation including CUDA GPU routines"
    if [[ "$2" == "circe" ]]; then
        echo "... for CIRCE cluster environment.";
        module purge
        module load compilers/gcc/4.9.4 # CUDA 7.5 not compatible with gcc > 5.0
        module load apps/cuda/7.5
        nvcc -x cu main.cpp -std=c++11 -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES -D__STRICT_ANSI__ -D CUDA -O3 -o ../mcmd
    elif [[ "$2" == "debug" ]]; then
        echo "... in debug mode (with errors/warnings)";
        nvcc -x cu main.cpp -std=c++11 -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES -D__STRICT_ANSI__ -D CUDA -G -g -O3 -o ../mcmd
    else
        nvcc -x cu main.cpp -std=c++11 -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES -D__STRICT_ANSI__ -D CUDA -O3 -o ../mcmd
    fi
########################################################################
elif [[ "$option" == "icpu" ]]; then
    # CPU compilation using Intel
    echo "Doing serial Intel (1 proc.) compilation for CPU"
    if [[ "$2" == "bridges" ]]; then
        echo "... for Bridges cluster environment";
        module purge
        module load icc/16.0.3
        icpc --std=c++11 -fast -unroll-aggressive -O3 -o ../mcmd main.cpp
    else
        icpc --std=c++11 -fast -unroll-aggressive -O3 -o ../mcmd main.cpp
    fi
#######################################################################
elif [[ "$option" == "mpi" ]]; then
    # MPI compilation
    echo "MPI is not implemented in MCMD as of now."
    #echo "Doing MPI compilation (for parallel implementation)"
    #mpic++ main.cpp -lm -o ../mcmd -I. -std=c++11 -D MPI -Ofast
#######################################################################
elif [[ "$option" == "omp" ]]; then 
    echo "Doing OpenMP compilation"
    if [[ "$2" == "linux" ]]; then
        echo "... for linux."
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -foptimize-sibling-calls -finline-limit=10000 -fexpensive-optimizations -flto -march=native -frename-registers -fopenmp -D OMP
    else 
        /usr/bin/llvm-g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -fopenmp
    fi
fi

echo "...done. Have a nice day."
