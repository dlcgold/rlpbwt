cmake_minimum_required(VERSION 3.20)
project(rlpbwt)

set(CMAKE_CXX_STANDARD 11)
add_compile_options(-Wall -Wextra -pedantic -O3)
add_executable(rlpbwt main.cpp lib/rlpbwt.cpp include/rlpbwt.h lib/rlpbwt_column.cpp include/rlpbwt_column.h lib/rlpbwt_rlrow.cpp include/rlpbwt_rlrow.h lib/utils.cpp include/utils.h include/exceptions.h lib/rlpbwt_match.cpp include/rlpbwt_match.h lib/rlpbwtc.cpp include/rlpbwtc.h lib/rlpbwtm.cpp include/rlpbwtm.h)
target_link_libraries(rlpbwt -lsdsl -ldivsufsort -ldivsufsort64 -lgtest)
#target_link_libraries(rlpbwt -lgtest)
