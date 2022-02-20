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

message("${shaped_slp_SOURCE_DIR}")
message("${FOLCA_SOURCE_DIR}")
message("${SUX_SOURCE_DIR}/function")
message("${SUX_SOURCE_DIR}/support")
message("${shaped_slp_SOURCE_DIR}/external/sux/sux")

option(BUILD_TESTS "Set OFF to not compile the tests" ON)
if (${BUILD_TESTS})
    enable_testing()
    add_executable(rlpbwt_test test/test.cpp
            lib/column.cpp include/column.h
            lib/utils.cpp include/utils.h
            include/exceptions.h
            lib/rlpbwt.cpp include/rlpbwt.h
            lib/match.cpp include/match.h
            lib/rlrow.cpp include/rlrow.h
            lib/birlpbwt.cpp include/birlpbwt.h
            lib/match_end.cpp include/match_end.h
            lib/rlpbwtbv.cpp include/rlpbwtbv.h
            lib/columnbv.cpp include/columnbv.h
            lib/column_thr.cpp include/column_thr.h
            lib/panel_ra.cpp include/panel_ra.h
            lib/slp_panel_ra.cpp include/slp_panel_ra.h
            include/rlpbwt_ra.h
            include/phi_support.h
            lib/ms.cpp include/ms.h
            lib/ms_matches.cpp include/ms_matches.h)

    #target_compile_options(rlpbwt_test PUBLIC "-std=c++17")
    #find_library(SDSL_lib lsdsl)
    #if (NOT SDSL_lib)
    #else ()
    #    target_link_libraries(rlpbwt_test -lsdsl -ldivsufsort -ldivsufsort64 -lgtest)
    #endif ()
    target_include_directories(rlpbwt_test PUBLIC
            "${shaped_slp_SOURCE_DIR}"
            "${FOLCA_SOURCE_DIR}"
            "${SUX_SOURCE_DIR}/function"
            "${SUX_SOURCE_DIR}/support")

    target_link_libraries(rlpbwt_test sdsl divsufsort divsufsort64 -lgtest)

    #target_link_libraries(rlpbwt_test -lsdsl -ldivsufsort -ldivsufsort64 -lgtest)
    add_test(NAME rlpbwt_test COMMAND rlpbwt_test)
else ()
    add_executable(rlpbwt build_rlpbwt_thr.cpp
            lib/column.cpp include/column.h
            lib/utils.cpp include/utils.h
            include/exceptions.h
            lib/rlpbwt.cpp include/rlpbwt.h
            lib/match.cpp include/match.h
            lib/rlrow.cpp include/rlrow.h
            lib/birlpbwt.cpp include/birlpbwt.h
            lib/match_end.cpp include/match_end.h
            lib/rlpbwtbv.cpp include/rlpbwtbv.h
            lib/columnbv.cpp include/columnbv.h
            lib/panel_ra.cpp include/panel_ra.h
            lib/slp_panel_ra.cpp include/slp_panel_ra.h
            lib/column_thr.cpp include/column_thr.h
            include/rlpbwt_ra.h
            include/phi_support.h
            lib/ms.cpp include/ms.h
            lib/ms_matches.cpp include/ms_matches.h)
    #target_compile_options(rlpbwt PUBLIC "-std=c++17")
    #find_library(SDSL_lib lsdsl)
    #if (NOT SDSL_lib)
    #target_link_libraries(rlpbwt sdsl divsufsort divsufsort64)
    #else ()
    #   target_link_libraries(rlpbwt -lsdsl -ldivsufsort -ldivsufsort64 -lgtest)
    #endif ()
    target_include_directories(rlpbwt PUBLIC
            "${shaped_slp_SOURCE_DIR}"
            "${FOLCA_SOURCE_DIR}"
            "${SUX_SOURCE_DIR}/function"
            "${SUX_SOURCE_DIR}/support"
            )
    target_link_libraries(rlpbwt sdsl divsufsort divsufsort64)
endif ()
