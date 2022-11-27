cmake_minimum_required(VERSION 3.2.0)
project(ppx VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 14)
add_executable(ppx_main main.cpp)

set(CMAKE_CXX_FLAGS_RELEASE "/MT")
set(CMAKE_CXX_FLAGS_DEBUG "/MTd")

include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.12.0
)
FetchContent_MakeAvailable(googletest)
add_executable(ppx_test test.cpp)
target_link_libraries(ppx_test gtest_main)

FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.7.1
)
FetchContent_MakeAvailable(benchmark)

# find_package(benchmark REQUIRED)
# add_executable(ppx_bench modern_robotics.cpp benchmark.cpp)
# target_include_directories(ppx_bench PRIVATE Eigen)
# target_link_libraries(ppx_bench benchmark::benchmark)
