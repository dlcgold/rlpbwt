cmake_minimum_required(VERSION 3.20)
project(rlpbwt)


include(FetchContent)

# External Project resources
# ------------------------------------------------------------------------------

# DIVSUFSORT
find_library(DIVSUFSORT_LIB divsufsort HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
if(NOT DIVSUFSORT_LIB)
    message(STATUS "divsufsort library not found. Building as an external content")
else()
    message(STATUS "divsufsort library found at ${DIVSUFSORT_LIB}.")
endif()

# DIVSUFSORT64
find_library(DIVSUFSORT64_LIB divsufsort64 HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
if(NOT DIVSUFSORT64_LIB)
    message(STATUS "divsufsort64 library not found. Building as an external content")
else()
    message(STATUS "divsufsort64 library found at ${DIVSUFSORT64_LIB}.")
endif()

# SDSL
find_library(SDSL_LIB sdsl HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
find_path(SDSL_SRC sdsl HINTS ${CMAKE_CURRENT_BINARY_DIR}/include PATHS ${CONDA_PREFIX}/include)
if(NOT SDSL_LIB)
    message(STATUS "sdsl library not found. Building as an external content")
    FetchContent_Declare(sdsl
            GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
            # GIT_TAG "v2.1.1"
            )

    if(NOT sdsl_POPULATED)
        FetchContent_Populate(sdsl)
        # Install in the current binary dir
        execute_process(COMMAND ./install.sh ${CMAKE_CURRENT_BINARY_DIR}
                RESULT_VARIABLE result
                WORKING_DIRECTORY ${sdsl_SOURCE_DIR}
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log)
        file(APPEND ${INSTALL_FILE} "MESSAGE(STATUS \"Installing sdsl.\")\n")
        file(APPEND ${INSTALL_FILE} "execute_process(COMMAND ./install.sh \${CMAKE_INSTALL_PREFIX}
                WORKING_DIRECTORY ${sdsl_SOURCE_DIR}
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log)\n")
        if(result)
            message(FATAL_ERROR "CMake step for sdsl-lite failed: ${result}")
        endif()
    endif()

    set(SDSL_SRC ${CMAKE_CURRENT_BINARY_DIR}/include)
    set(DIVSUFSORT_SRC ${CMAKE_CURRENT_BINARY_DIR}/include)
    set(SDSL_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libsdsl.a)
    set(DIVSUFSORT_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libdivsufsort.a)
    set(DIVSUFSORT64_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libdivsufsort64.a)

else()
    message(STATUS "sdsl library found at ${SDSL_LIB}.")
    message(STATUS "sdsl sources found at ${SDSL_SRC}.")
endif()

add_library(sdsl STATIC IMPORTED GLOBAL)
set_property(TARGET sdsl PROPERTY IMPORTED_LOCATION ${SDSL_LIB})
target_include_directories(sdsl INTERFACE ${SDSL_SRC})

add_library(divsufsort STATIC IMPORTED GLOBAL)
set_property(TARGET divsufsort PROPERTY IMPORTED_LOCATION ${DIVSUFSORT_LIB})

add_library(divsufsort64 STATIC IMPORTED GLOBAL)
set_property(TARGET divsufsort64 PROPERTY IMPORTED_LOCATION ${DIVSUFSORT64_LIB})

set(CMAKE_CXX_STANDARD 14)
add_compile_options(-Wall -Wextra -pedantic -O3)
add_executable(rlpbwt build_rlpbwt_thr.cpp lib/column.cpp include/column.h lib/utils.cpp include/utils.h include/exceptions.h lib/rlpbwt.cpp include/rlpbwt.h lib/match.cpp include/match.h lib/rlrow.cpp include/rlrow.h lib/birlpbwt.cpp include/birlpbwt.h lib/match_end.cpp include/match_end.h lib/rlpbwtbv.cpp include/rlpbwtbv.h lib/columnbv.cpp include/columnbv.h lib/rlpbwt_thr.cpp include/rlpbwt_thr.h lib/column_thr.cpp include/column_thr.h lib/panel_ra.cpp include/panel_ra.h)
#add_executable(rlpbwt main.cpp lib/column.cpp include/column.h lib/utils.cpp include/utils.h include/exceptions.h lib/rlpbwt.cpp include/rlpbwt.h lib/match.cpp include/match.h lib/rlrow.cpp include/rlrow.h lib/birlpbwt.cpp include/birlpbwt.h lib/match_end.cpp include/match_end.h lib/rlpbwtbv.cpp include/rlpbwtbv.h lib/columnbv.cpp include/columnbv.h lib/rlpbwt_thr.cpp include/rlpbwt_thr.h lib/column_thr.cpp include/column_thr.h lib/panel_ra.cpp include/panel_ra.h)
target_link_libraries(rlpbwt -lsdsl -ldivsufsort -ldivsufsort64 -lgtest -ltcmalloc -lboost_serialization -lboost_system)
