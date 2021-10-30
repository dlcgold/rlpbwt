cmake_minimum_required(VERSION 3.20)
project(rlpbwt)

set(CMAKE_CXX_STANDARD 11)
add_compile_options(-Wall -Wextra -pedantic)
add_executable(rlpbwt main.cpp lib/rlpbwt.cpp include/rlpbwt.h lib/pbwt_column.cpp include/pbwt_column.h lib/pbwt_rlrow.cpp include/pbwt_rlrow.h lib/utils.cpp include/utils.h include/exceptions.h)
#target_link_libraries(rlpbwt -lsdsl -ldivsufsort -ldivsufsort64)
