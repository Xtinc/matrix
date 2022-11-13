cmake_minimum_required(VERSION 3.0.0)
project(matrix VERSION 0.1.0)

option(ENABLE_BENCHMARK OFF)

if(ENABLE_BENCHMARK)
    find_package(benchmark REQUIRED)
    add_executable(matrix_bench benchmark.cpp)
    target_include_directories(matrix_bench PRIVATE Eigen)
endif()

add_executable(matrix_test main.cpp)
