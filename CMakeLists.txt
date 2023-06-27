cmake_minimum_required(VERSION 3.11)
project(ppx VERSION 0.1.0)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

if(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W3" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_library(ppx INTERFACE)
target_include_directories(ppx INTERFACE include)

add_subdirectory(demo)

option(ENABLE_UNIT_TESTS "Enable unit tests" ON)
option(ENABLE_LCOV "Enable code coverage tests" OFF)
message(STATUS "---------------------------")
message(STATUS "Current : ${PROJECT_NAME}")
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
# add_executable(ppx_bench misc/benchmark.cpp)
# target_include_directories(ppx_bench PRIVATE include Eigen)
# target_link_libraries(ppx_bench benchmark::benchmark)