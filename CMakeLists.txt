cmake_minimum_required(VERSION 3.20)
project(rlpbwt)

set(CMAKE_CXX_STANDARD 11)
add_compile_options(-Wall -Wextra -pedantic -O3)
#add_executable(rlpbwt main.cpp lib/rlpbwt.cpp include/rlpbwt.h lib/column.cpp include/column.h lib/rlpbwt_rlrow.cpp include/rlpbwt_rlrow.h lib/utils.cpp include/utils.h include/exceptions.h lib/rlpbwt_match.cpp include/rlpbwt_match.h lib/rlpbwt.cpp include/rlpbwt.h lib/match.cpp include/match.h lib/rlrow.cpp include/rlrow.h)
add_executable(rlpbwt main.cpp lib/column.cpp include/column.h lib/utils.cpp include/utils.h include/exceptions.h lib/rlpbwt.cpp include/rlpbwt.h lib/match.cpp include/match.h lib/rlrow.cpp include/rlrow.h lib/birlpbwt.cpp include/birlpbwt.h)
target_link_libraries(rlpbwt -lsdsl -ldivsufsort -ldivsufsort64 -lgtest -ltcmalloc)
#target_link_libraries(rlpbwt -lgtest)
