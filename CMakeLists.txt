cmake_minimum_required(VERSION 3.8)
project(mcmd)
set(CMAKE_CXX_FLAGS "-g3")

set(CMAKE_CXX_STANDARD 11)

set(INCLUDE src/)

set(SOURCE_FILES src/main.cpp )


add_executable(mcmd ${SOURCE_FILES} ${INCLUDE})