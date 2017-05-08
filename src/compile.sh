#!/bin/bash

#########################
# HOW TO USE THIS SCRIPT
# bash compile.sh [type] [HPC cluster if any]
#
# e.g.
#
# bash compile.sh cpu           (for a single computer on Mac or Linux)
# bash compile.sh cpu linux     (optimized for linux)
# bash compile.sh cpu errors    (same but with a lot of errors/warnings)
# bash compile.sh cpu circe     (for CIRCE)
# bash compile.sh cpu bridges   (for bridges)
# bash compile.sh icpu bridges  (for bridges Intel compilation (seems slower than gcc))
# bash compile.sh gpu           (for a single computer with GPU functions with NVIDIA CUDA installed).
# bash compile.sh gpu circe     (GPU on circe, still not working)
##########################

echo "This should take like 10 sec."

DEFAULT="cpu"

if [ $# -eq 0 ]; then # IF NO ARGUMENT GIVEN
option=$DEFAULT   # MANUAL OPTION (only reads if no argument given)
elif [ $1 == "cpu" ]; then
option="cpu"
elif [ $1 == "gpu" ]; then
option="gpu";
fi


if [[ "$option" == "cpu" ]]; then
    # THIS IS FOR SERIAL COMPILATION (1 CPU ONLY, NO GPU)
    echo "Doing serial GCC (1 processor) compilation for CPU"
    if [[ "$2" == "circe" ]]; then
        module purge
        module load compilers/gcc/6.2.0
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast;
    elif [[ "$2" == "bridges" ]]; then
        module purge
        module load gcc/6.3.0
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast;
    elif [[ "$2" == "errors" ]]; then
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -Werror -Wall;
    elif [[ "$2" == "linux" ]]; then
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -foptimize-sibling-calls -finline-limit=10000 -fexpensive-optimizations -flto -march=native -frename-registers 
    else
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast;
    fi

elif [[ "$option" == "gpu" ]]; then
    # GPU compilation enabled
    echo "Doing serial GCC (1 processor) compilation including CUDA GPU routines"
    if [[ "$2" == "circe" ]]; then
        module purge
        module load compilers/gcc/6.2.0
        module load apps/cuda/7.5
        nvcc -x cu main.cpp -std=c++11 -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES -D__STRICT_ANSI__ -D CUDA -O3 -o ../mcmd
    else
        nvcc -x cu main.cpp -std=c++11 -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES -D__STRICT_ANSI__ -D CUDA -O3 -o ../mcmd
    fi
elif [[ "$option" == "icpu" ]]; then
    # CPU compilation using Intel
    echo "Doing serial Intel (1 proc.) compilation for CPU"
    if [[ "$2" == "bridges" ]]; then
        module purge
        module load icc/16.0.3
        icpc --std=c++11 -fast -unroll-aggressive -O3 -o ../mcmd main.cpp
    else
        icpc --std=c++11 -fast -unroll-aggressive -O3 -o ../mcmd main.cpp
    fi
fi
