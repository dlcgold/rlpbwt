cmake_minimum_required(VERSION 3.16)
project(RLPBWT)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -pedantic -O3)

find_package(OpenMP)

# Configure thirdparty
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_INCLUDEDIR "include")

add_subdirectory(thirdparty)

FetchContent_GetProperties(shaped_slp)
set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)

option(BUILD_TESTS "Set OFF to not compile the tests" ON)

if (${BUILD_TESTS})
    find_package(benchmark REQUIRED)
    enable_testing()
    add_executable(rlpbwt_test test/test.cpp
            lib/utils.cpp include/utils.h
            include/exceptions.h
            lib/rlpbwt_bv.cpp include/rlpbwt_bv.h
            lib/column_bv.cpp include/column_bv.h
            lib/column_ms.cpp include/column_ms.h
            lib/panel_ra.cpp include/panel_ra.h
            lib/slp_panel_ra.cpp include/slp_panel_ra.h
            include/rlpbwt_ms.h
            include/phi_support.h
            lib/ms.cpp include/ms.h
            lib/ms_matches.cpp include/ms_matches.h
            lib/column_naive.cpp include/column_naive.h
            lib/rlpbwt_naive.cpp include/rlpbwt_naive.h
            lib/matches_naive.cpp include/matches_naive.h
            test/benchmark_test.h)

    target_include_directories(rlpbwt_test PUBLIC
            "${shaped_slp_SOURCE_DIR}"
            "${FOLCA_SOURCE_DIR}"
            "${SUX_SOURCE_DIR}/function"
            "${SUX_SOURCE_DIR}/support")

    target_link_libraries(rlpbwt_test sdsl divsufsort divsufsort64 -lgtest
            benchmark::benchmark)

    add_test(NAME rlpbwt_test COMMAND rlpbwt_test)
else ()
    add_executable(rlpbwt rlpbwt_main.cpp
            lib/utils.cpp include/utils.h
            include/exceptions.h
            lib/rlpbwt_bv.cpp include/rlpbwt_bv.h
            lib/column_bv.cpp include/column_bv.h
            lib/panel_ra.cpp include/panel_ra.h
            lib/slp_panel_ra.cpp include/slp_panel_ra.h
            lib/column_ms.cpp include/column_ms.h
            include/rlpbwt_ms.h
            include/phi_support.h
            lib/ms.cpp include/ms.h
            lib/ms_matches.cpp include/ms_matches.h
            lib/column_naive.cpp include/column_naive.h
            lib/rlpbwt_naive.cpp include/rlpbwt_naive.h
            lib/matches_naive.cpp include/matches_naive.h)

    target_include_directories(rlpbwt PUBLIC
            "${shaped_slp_SOURCE_DIR}"
            "${FOLCA_SOURCE_DIR}"
            "${SUX_SOURCE_DIR}/function"
            "${SUX_SOURCE_DIR}/support"
            )
    target_link_libraries(rlpbwt sdsl divsufsort divsufsort64)
endif ()
