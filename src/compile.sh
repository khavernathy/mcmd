#!/bin/bash


g++ main.cpp -lm -o ../t -I. -std=c++11; 
cp ../t ../md_test/; 
cp ../t ../npt_test/; 
cp ../t ../uvt_test;
cp ../t ../nve_test;
cp ../t ../nopbc_test
