cmake_minimum_required(VERSION 3.0.0)
project(matrix VERSION 0.1.0)

find_package(benchmark REQUIRED)
add_executable(matrix_bench benchmark.cpp)
target_include_directories(matrix_bench PRIVATE Eigen)
target_link_libraries(matrix_bench benchmark::benchmark)
add_executable(matrix_test main.cpp)
