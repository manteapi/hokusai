cmake_minimum_required(VERSION 2.6)
PROJECT(HOKUSAI)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-std=c++0x -fopenmp" ${CMAKE_CXX_CFLAGS})

set(HOKUSAI_SOURCE
    ./src/particle.cpp
    ./src/system.cpp
    ./src/utils.cpp
    ./src/HSL2RGB.cpp
    ./src/write_bmp.cpp
    ./src/gridUtility.cpp
    ./src/io.cpp
    ./src/rasterizer.cpp
    ./src/triMesh.cpp
    ./src/sampler.cpp
)

set(HOKUSAI_HEADER
    ./include/hokusai/particle.hpp
    ./include/hokusai/system.hpp
    ./include/hokusai/utils.hpp
    ./include/hokusai/HSL2RGB.hpp
    ./include/hokusai/Vec.hpp
    ./include/hokusai/write_bmp.hpp
    ./include/hokusai/gridUtility.hpp
    ./include/hokusai/io.hpp
    ./include/hokusai/rasterizer.hpp
    ./include/hokusai/triMesh.hpp
    ./include/hokusai/strutil.hpp
    ./include/hokusai/sampler.hpp
)

include_directories( include )

add_library(HOKUSAI STATIC
    ${HOKUSAI_SOURCE}
    ${HOKUSAI_HEADER}
)
