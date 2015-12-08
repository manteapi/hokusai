cmake_minimum_required(VERSION 2.6)
project(HOKUSAI) #Make HOKUSAI_SOURCE_DIR, HOKUSAI_BIN_DIR and HOKUSAI available

#set(CMAKE_BUILD_TYPE Release)
##set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    option(HOKUSAI_OPENMP "Use OpenMP" ON)
else()
    option(HOKUSAI_OPENMP "Use OpenMP" OFF)
endif()

if(HOKUSAI_OPENMP)
    add_definitions( -DHOKUSAI_USING_OPENMP )
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

set(HOKUSAI_SOURCE
    ./src/particle.cpp
    ./src/system.cpp
    ./src/utils.cpp
    ./src/gridUtility.cpp
    ./src/marchingSquare.cpp
    ./src/implicitFunction.cpp
    ./src/io.cpp
    ./src/particleSource.cpp
    ./src/particleSink.cpp
)

set(HOKUSAI_HEADER
    ./include/hokusai/utility.hpp
    ./include/hokusai/particle.hpp
    ./include/hokusai/system.hpp
    ./include/hokusai/utils.hpp
    ./include/hokusai/gridUtility.hpp
    ./include/hokusai/marchingSquare.hpp
    ./include/hokusai/implicitFunction.hpp
    ./include/hokusai/io.hpp
    ./include/hokusai/particleSource.hpp
    ./include/hokusai/particleSink.hpp
)

set(MAGNET_INCLUDE_DIRS "" CACHE PATH "path to magnet headers")
include_directories(${MAGNET_INCLUDE_DIRS})

set(EIGEN_INCLUDE_DIRS "" CACHE PATH "path to eigen headers")
include_directories(${EIGEN_INCLUDE_DIRS})

include_directories( include )

find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        file(MAKE_DIRECTORY doc) 
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif(DOXYGEN_FOUND)

set(EXTRA_CLEANING
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/html
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/latex
    )
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${EXTRA_CLEANING}")

add_library(HOKUSAI STATIC ${HOKUSAI_SOURCE} ${HOKUSAI_HEADER})

add_executable(benchmark_sink benchmark/benchmark_sink.cpp ${HOKUSAI_SOURCE} ${HOKUSAI_HEADER})
add_executable(benchmark_source benchmark/benchmark_source.cpp ${HOKUSAI_SOURCE} ${HOKUSAI_HEADER})
add_executable(benchmark_damBreak benchmark/benchmark_damBreak.cpp ${HOKUSAI_SOURCE} ${HOKUSAI_HEADER})
add_executable(benchmark_splash benchmark/benchmark_splash.cpp ${HOKUSAI_SOURCE} ${HOKUSAI_HEADER})
add_executable(benchmark_wave benchmark/benchmark_wave.cpp ${HOKUSAI_SOURCE} ${HOKUSAI_HEADER})
