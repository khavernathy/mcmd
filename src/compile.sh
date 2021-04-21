#!/bin/bash

#########################
# HOW TO USE THIS SCRIPT
# bash compile.sh [type] [HPC cluster if any]
#
# e.g.
#
# bash compile.sh               Defaults to cpu compilation
# bash compile.sh cpu           (same result; for a single computer on Mac or Linux or RaspPi)
# bash compile.sh cpu linux     (optimized for linux (ONLY))
# bash compile.sh cpu windows   (for use on Windows ONLY -- you must have gcc installed, e.g. through Cygwin)
	# note -- to install on windows through the command line,
	# run the following commands:
	# 	"C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\Tools\vsdevcmd"
	# 	cl /EHsc main.cpp /D "WINDOWS" /Fe..\mcmd
	# this will write a windows executable called mcmd.exe in the main directory
# bash compile.sh cpu debug     (same but with errors/warnings)
# bash compile.sh cpu circe     (for CIRCE (HPC at USF))
# bash compile.sh cpu bridges   (for bridges (HPC at UPitt)

# bash compile.sh icpu bridges  (for bridges Intel compilation (seems slower than gcc))

# bash compile.sh cmake         (compile serial version using cmake. Not recommended for best performance.)

# bash compile.sh gpu           (for a single computer with CUDA compatible GPU installed).
# bash compile.sh gpu debug     (for GPU compilation with errors)
# bash compile.sh gpu circe     (CUDA GPU on circe HPC)
# bash compile.sh gpu windows   (CUDA GPU on Windows machine)

# bash compile.sh omp           (for OpenMP implementation)
# bash compile.sh omp linux     (same, optimized for Linux)

# bash compile.sh mpi           (for MPI implementation, not yet functional)


##########################


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
elif [ $1 == "cmake" ]; then
    option="cmake"
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
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -foptimize-sibling-calls -finline-limit=10000 -fexpensive-optimizations -flto  -frename-registers
    elif [[ "$2" == "bridges" ]]; then
        echo "... for Bridges cluster environment.";
        module purge
        module load gcc/6.3.0
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -foptimize-sibling-calls -finline-limit=10000 -fexpensive-optimizations -flto -march=native -frename-registers
    elif [[ "$2" == "debug" ]]; then
        echo "... in debug mode (with errors/warnings";
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -Werror -Wall -ggdb;
    elif [[ "$2" == "linux" ]]; then
        echo "... optimized for linux."
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -foptimize-sibling-calls -finline-limit=10000 -fexpensive-optimizations -flto -march=native -frename-registers 
    elif [[ "$2" == "O3" ]]; then
        echo "... optimized at -O3."
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -O3
    elif [[ "$2" == "windows" ]]; then
        echo "... for Windows OS.";
	echo "Note: to install in the Windows command line, use the following commands:"
	echo "C:\\Program Files (x86)\\Microsoft Visual Studio\\2017\\Community\\Common7\\Tools\\vsdevcmd"
	echo "cl /EHsc main.cpp /D \"WINDOWS\" /Fe..\\mcmd"
        echo ""
	echo "=================================================="
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
    elif [[ "$2" == "windows" ]]; then
	echo "... for Windows OS."
	nvcc -x cu main.cpp -std=c++11 -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES -D__STRICT_ANSI__ -D CUDA -O3 -D WINDOWS -o ../mcmd
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
        icpc --std=c++11 -fast -unroll-aggressive -Ofast -o ../mcmd main.cpp
    else
        icpc --std=c++11 -fast -unroll-aggressive -Ofast -o ../mcmd main.cpp
    fi
#######################################################################
elif [[ "$option" == "cmake" ]]; then
    # cmake serial compilation
    echo "Using cmake to compile MCMD."
    cd ..
    cmake . -DCMAKE_BUILD_TYPE=Release
    make -j12
    cd src
#######################################################################
elif [[ "$option" == "mpi" ]]; then
    # MPI compilation
    echo "MPI is not implemented in MCMD as of now."
    #echo "Doing MPI compilation (for parallel implementation)"
    #mpic++ main.cpp -lm -o ../mcmd -I. -std=c++11 -D MPI -Ofast
#######################################################################
elif [[ "$option" == "omp" ]]; then 
    # OpenMP compilation
    #echo "OpenMP is not implemented in MCMD as of now."
    echo "Doing OpenMP compilation"
    if [[ "$2" == "linux" ]]; then
        echo "... for linux."
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -foptimize-sibling-calls -finline-limit=10000 -fexpensive-optimizations -flto -march=native -frename-registers -fopenmp -D OMP
    elif [[ "$2" == "circe" ]]; then
        echo "... for CIRCE cluster environment.";
        module purge
        module load compilers/gcc/6.2.0
        g++ main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -foptimize-sibling-calls -finline-limit=10000 -fexpensive-optimizations -flto -frename-registers -fopenmp -D OMP
    else 
        /usr/local/bin/g++-4.9 main.cpp -lm -o ../mcmd -I. -std=c++11 -Ofast -fopenmp -D OMP
    fi
fi

echo ""
echo ""
echo "... done. As long as you see no compilation errors above, you're good to go."
echo "... If you did, check the contents of compile.sh for examples"
echo "... or contact us by making an Issue on GitHub: https://github.com/khavernathy/mcmd/issues"
echo "... Stay awesome and have a nice day. :-)"
