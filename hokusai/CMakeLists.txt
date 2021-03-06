cmake_minimum_required(VERSION 2.6)
project(HOKUSAI) #Make HOKUSAI_SOURCE_DIR, HOKUSAI_BIN_DIR and HOKUSAI available

set(CMAKE_BUILD_TYPE Release)
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
    ./src/kernel.cpp
    ./src/system.cpp
    ./src/particle.cpp
    ./src/boundary.cpp
    ./src/utils.cpp
    ./src/HSL2RGB.cpp
    ./src/write_bmp.cpp
    ./src/gridUtility.cpp
    ./src/io.cpp
    ./src/rasterizer.cpp
    ./src/triMesh.cpp
    ./src/sampler.cpp
    ./src/particleSource.cpp
    ./src/particleSink.cpp
    ./src/fluidParams.cpp
    ./src/boundaryParams.cpp
    ./src/solverParams.cpp
)

set(HOKUSAI_HEADER
    ./include/hokusai/kernel.hpp
    ./include/hokusai/system.hpp
    ./include/hokusai/common.hpp
    ./include/hokusai/particle.hpp
    ./include/hokusai/boundary.hpp
    ./include/hokusai/utils.hpp
    ./include/hokusai/HSL2RGB.hpp
    ./include/hokusai/write_bmp.hpp
    ./include/hokusai/gridUtility.hpp
    ./include/hokusai/io.hpp
    ./include/hokusai/rasterizer.hpp
    ./include/hokusai/triMesh.hpp
    ./include/hokusai/strutil.hpp
    ./include/hokusai/sampler.hpp
    ./include/hokusai/particleSource.hpp
    ./include/hokusai/particleSink.hpp
    ./include/hokusai/fluidParams.hpp
    ./include/hokusai/boundaryParams.hpp
    ./include/hokusai/solverParams.hpp
)

set(ALJABR_INCLUDE_DIRS "" CACHE PATH "path to aljabr headers")
include_directories(${ALJABR_INCLUDE_DIRS})

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

add_library(HOKUSAI STATIC
    ${HOKUSAI_SOURCE}
    ${HOKUSAI_HEADER}
)
