cmake_minimum_required(VERSION 3.20)
project(cspike)

set(CMAKE_CXX_STANDARD 11)

add_compile_options(-std=c++11 -O3 -m64)
add_executable(cspike_std_mem main.cpp utils/MurmurHash3.cpp spike_sketch_extend.h  impl/ss_query3.cpp)






