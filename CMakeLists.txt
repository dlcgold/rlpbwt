cmake_minimum_required(VERSION 3.20)
project(rlpbwt)

set(CMAKE_CXX_STANDARD 14)
add_compile_options(-Wall -Wextra -pedantic -O3)

# Configure thirdparty
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_INCLUDEDIR "include")

add_subdirectory(thirdparty)
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
            lib/rlpbwt_thr.cpp include/rlpbwt_thr.h
            lib/column_thr.cpp include/column_thr.h
            lib/panel_ra.cpp include/panel_ra.h)

    find_library(SDSL_lib lsdsl)
    if (NOT SDSL_lib)
        target_link_libraries(rlpbwt_test sdsl divsufsort divsufsort64 -lgtest)
    else ()
        target_link_libraries(rlpbwt_test -lsdsl -ldivsufsort -ldivsufsort64 -lgtest)
    endif ()
    target_link_libraries(rlpbwt_test -lsdsl -ldivsufsort -ldivsufsort64 -lgtest)
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
            lib/rlpbwt_thr.cpp include/rlpbwt_thr.h
            lib/column_thr.cpp include/column_thr.h
            lib/panel_ra.cpp include/panel_ra.h)

    if (NOT SDSL_lib)
        target_link_libraries(rlpbwt sdsl divsufsort divsufsort64 -lgtest)
    else ()
        target_link_libraries(rlpbwt -lsdsl -ldivsufsort -ldivsufsort64 -lgtest)
    endif ()
endif ()
