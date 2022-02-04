
cmake_minimum_required(VERSION 3.10)
project(rlpbwt)

set(CMAKE_CXX_STANDARD 14)
add_compile_options(-Wall -Wextra -pedantic -O3)

# Configure thirdparty
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_INCLUDEDIR "include") # This is an hack because include(GUIInstallDirs) doesn't work

add_subdirectory(thirdparty)

add_executable(rlpbwt_build build_rlpbwt_thr.cpp lib/utils.cpp include/utils.h include/exceptions.h lib/rlpbwt_thr.cpp include/rlpbwt_thr.h lib/column_thr.cpp include/column_thr.h lib/panel_ra.cpp include/panel_ra.h)
target_link_libraries(rlpbwt_build sdsl divsufsort divsufsort64)

