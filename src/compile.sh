#!/bin/bash

#########################
# HOW TO USE THIS SCRIPT
# Make sure you have the needed compiler loaded first. (gcc 6.3.0 or 6.2.0 etc.)
#
#
# bash compile.sh cpu
# ---- or ----
# bash compile.sh gpu    (for GPU you will need to be on a GPU cluster environment with CUDA loaded).
#
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
    g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast; 
elif [[ "$option" == "gpu" ]]; then
    # GPU compilation enabled
    echo "Doing serial GCC (1 processor) compilation including CUDA GPU routines"
    nvcc -x cu main.cpp -std=c++11 -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES -D__STRICT_ANSI__ -D CUDA -O3 -o ../mcmd
elif [[ "$option" == "icpu" ]]; then
    # CPU compilation using Intel
    echo "Doing serial Intel (1 proc.) compilation for CPU"
    icpc --std=c++11 -fast -unroll-aggressive -o ../mcmd main.cpp

fi
