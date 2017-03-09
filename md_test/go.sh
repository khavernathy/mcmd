#!/bin/bash 

cd ../src/
g++ main.cpp -lm -o ../t -I. -std=c++11
cd ../
cp t md_test
cd md_test
./t *inp | tee runlog
