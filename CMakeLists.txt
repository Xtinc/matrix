cmake_minimum_required(VERSION 3.11)
project(ppx VERSION 0.1.0)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories(include)

add_executable(ppx_demo demo/ppxlog.cpp demo/main.cpp)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_link_libraries(ppx_demo pthread)
endif()
if(MSVC)
  # target_compile_options(ppx_main PRIVATE /W4)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
else()
  target_compile_options(ppx_demo PRIVATE -Wall -Wextra -pedantic)
endif()


option(ENABLE_UNIT_TESTS "Enable unit tests" ON)
message(STATUS "---------------------------")
message(STATUS "Current Project: ${PROJECT_NAME}")
option(ENABLE_UNIT_TESTS "Enable unit tests" ON)
option(ENABLE_LCOV "Enable code coverage tests" OFF)
message(STATUS "Enable gtesting: ${ENABLE_UNIT_TESTS}")
message(STATUS "Enable coverage: ${ENABLE_LCOV}")
if(ENABLE_UNIT_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

# set(BENCHMARK_DOWNLOAD_DEPENDENCIES ON)
# include(FetchContent)

# FetchContent_Declare(
#     benchmark
#     GIT_REPOSITORY https://github.com/google/benchmark.git
#     GIT_TAG v1.7.1
# )
# FetchContent_MakeAvailable(benchmark)

# find_package(benchmark REQUIRED)
# add_executable(ppx_bench benchmark.cpp)
# target_include_directories(ppx_bench PRIVATE include)
# target_link_libraries(ppx_bench benchmark::benchmark)