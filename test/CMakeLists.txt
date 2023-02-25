set(SUBPRJ "test")
message(STATUS "---------------------------")
message(STATUS "Current : ${SUBPRJ}")
if(ENABLE_UNIT_TESTS)
message(STATUS "Enable Gtest in ${SUBPRJ}")
    message(STATUS "-- Gtest Configure")
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        release-1.12.0
    )
    if(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    endif()
    FetchContent_MakeAvailable(googletest)
    include(GoogleTest)
    add_executable(matrix_test matrix_test.cpp)
    target_link_libraries(matrix_test PRIVATE ppx gtest_main)
    # add_executable(bag_test bag_test.cpp)
    # target_link_libraries(bag_test PRIVATE BAGREC gtest_main)
    # add_executable(cpr_test cpr_test.cpp)
    # target_link_libraries(cpr_test PRIVATE CBOR gtest_main)
    gtest_discover_tests(matrix_test
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/build/test)
    # gtest_discover_tests( bag_test
    # WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/build/test)
    # gtest_discover_tests(cpr_test
    # WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/build/test)
    if(ENABLE_LCOV)
        message(STATUS "Enable Lcov in ${SUBPRJ}")
        message(STATUS "-- lcov Configure")
        find_program(LCOV_BIN lcov)
        if(NOT LCOV_BIN MATCHES "lcov$")
            message(FATAL_ERROR "lcov required, but not found!")
        endif()
        if(CMAKE_BUILD_TYPE MATCHES Debug)
            include(CodeCoverage.cmake)
            APPEND_COVERAGE_COMPILER_FLAGS()
            setup_target_for_coverage_lcov(
            NAME coverage
            EXECUTABLE ctest test
            BASE_DIRECTORY "${PROJECT_SOURCE_DIR}/build"
            EXCLUDE "${PROJECT_SOURCE_DIR}/build/_deps/*" "/usr/*"
            DEPENDENCIES matrix_test)
        endif()
    endif()
endif()