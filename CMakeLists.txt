cmake_minimum_required(VERSION 3.11)
project(ppx VERSION 1.2)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_AVX "Enable AVX intrinsic" OFF)
option(ENABLE_UNIT_TESTS "Enable unit tests" ON)
option(ENABLE_BENCHMARK "Enable benchmark" OFF)
option(ENABLE_LCOV "Enable code coverage tests" OFF)

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

  if(ENABLE_AVX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

  if(ENABLE_AVX AND COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mfma")
  endif()
endif()

add_library(ppx INTERFACE)
target_include_directories(ppx INTERFACE include)
if(ENABLE_AVX)
target_compile_definitions(ppx INTERFACE PPX_USE_AVX)
endif()
add_subdirectory(log)
add_subdirectory(demo)


if(ENABLE_UNIT_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

if(ENABLE_BENCHMARK)
  set(BENCHMARK_DOWNLOAD_DEPENDENCIES ON)
  include(FetchContent)
  FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.7.1
  )
  FetchContent_MakeAvailable(benchmark)
  add_subdirectory(benchmark)
endif()

message(STATUS "---------------------------")
message(STATUS "Current : ${PROJECT_NAME}")
message(STATUS "Enable AVX-INS  : ${ENABLE_AVX}")
message(STATUS "Enable gtesting : ${ENABLE_UNIT_TESTS}")
message(STATUS "Enable benchmark: ${ENABLE_BENCHMARK}")
message(STATUS "Enable coverage : ${ENABLE_LCOV}")
message(STATUS "---------------------------")