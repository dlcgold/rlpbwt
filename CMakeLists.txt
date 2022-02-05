cmake_minimum_required(VERSION 3.20)
project(rlpbwt)

set(CMAKE_CXX_STANDARD 14)
add_compile_options(-Wall -Wextra -pedantic -O3)

#add_executable(rlpbwt build_rlpbwt_thr.cpp lib/column.cpp include/column.h lib/utils.cpp include/utils.h include/exceptions.h lib/rlpbwt.cpp include/rlpbwt.h lib/match.cpp include/match.h lib/rlrow.cpp include/rlrow.h lib/birlpbwt.cpp include/birlpbwt.h lib/match_end.cpp include/match_end.h lib/rlpbwtbv.cpp include/rlpbwtbv.h lib/columnbv.cpp include/columnbv.h lib/rlpbwt_thr.cpp include/rlpbwt_thr.h lib/column_thr.cpp include/column_thr.h lib/panel_ra.cpp include/panel_ra.h)
add_executable(rlpbwt main.cpp lib/column.cpp include/column.h lib/utils.cpp include/utils.h include/exceptions.h lib/rlpbwt.cpp include/rlpbwt.h lib/match.cpp include/match.h lib/rlrow.cpp include/rlrow.h lib/birlpbwt.cpp include/birlpbwt.h lib/match_end.cpp include/match_end.h lib/rlpbwtbv.cpp include/rlpbwtbv.h lib/columnbv.cpp include/columnbv.h lib/rlpbwt_thr.cpp include/rlpbwt_thr.h lib/column_thr.cpp include/column_thr.h lib/panel_ra.cpp include/panel_ra.h)
target_link_libraries(rlpbwt -lsdsl -ldivsufsort -ldivsufsort64 -lgtest)
