set(SUBPRJ "test")
message(STATUS "---------------------------")
message(STATUS "Current : ${SUBPRJ}")
if (ENABLE_UNIT_TESTS)
    message(STATUS "Enable Gtest in ${SUBPRJ}")
    message(STATUS "-- Gtest Configure")
    include(FetchContent)
    FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG release-1.12.0
    )
    if (CMAKE_CXX_COMPILER_ID MATCHES MSVC)
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    endif ()
    FetchContent_MakeAvailable(googletest)
    include(GoogleTest)
    add_executable(matrix_test matrix_test.cpp)
    target_link_libraries(matrix_test PRIVATE ppx gtest_main)
    add_executable(lineqn_test solver_test.cpp)
    target_link_libraries(lineqn_test PRIVATE ppx gtest_main)
    add_executable(opt_test opt_test.cpp)
    target_link_libraries(opt_test PRIVATE ppx gtest_main)
    add_executable(lie_test lie_test.cpp)
    target_link_libraries(lie_test PRIVATE ppx gtest_main)
    add_executable(statistic_test statistic_test.cpp)
    target_link_libraries(statistic_test PRIVATE ppx gtest_main)
    gtest_discover_tests(matrix_test
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/build/test)
    gtest_discover_tests(lineqn_test
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/build/test)
    gtest_discover_tests(opt_test
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/build/test)
    gtest_discover_tests(lie_test
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/build/test)
    gtest_discover_tests(statistic_test
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/build/test)
    if (ENABLE_LCOV)
        message(STATUS "Enable Lcov in ${SUBPRJ}")
        message(STATUS "-- lcov Configure")
        find_program(LCOV_BIN lcov)
        if (NOT LCOV_BIN MATCHES "lcov$")
            message(FATAL_ERROR "lcov required, but not found!")
        endif ()
        if (CMAKE_BUILD_TYPE MATCHES Debug)
            include(CodeCoverage.cmake)
            APPEND_COVERAGE_COMPILER_FLAGS()
            setup_target_for_coverage_lcov(
                    NAME coverage
                    EXECUTABLE ctest test
                    BASE_DIRECTORY "${PROJECT_SOURCE_DIR}/build"
                    EXCLUDE "${PROJECT_SOURCE_DIR}/build/_deps/*" "/usr/*"
                    DEPENDENCIES matrix_test lineqn_test opt_test lie_test statistic_test)
        endif ()
    endif ()
endif ()