cmake_minimum_required(VERSION 2.6)
project(hokusai_examples)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-dev")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_subdirectory(../hokusai hokusai/build)
include_directories(../hokusai/include)

find_package(Boost COMPONENTS timer system REQUIRED)
include_directories( ${BOOST_INCLUDE_DIRS} )

set(ALJABR_INCLUDE_DIRS "" CACHE PATH "path to alajabr headers")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/eigen-3.3.4)

include_directories(${ALJABR_INCLUDE_DIRS})

if(HOKUSAI_OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

set(
    ALL_LIBS
    ${Boost_LIBRARIES}
    HOKUSAI 
    )

add_executable(mesher mesher.cpp rchandra/CIsoSurface.cpp rchandra/CIsoSurface.h rchandra/Vectors.cpp rchandra/Vectors.h)
target_link_libraries(mesher ${ALL_LIBS})

add_executable(particleSource particleSource.cpp)
target_link_libraries(particleSource ${ALL_LIBS})

add_executable(boundaryMesh boundaryMesh.cpp)
target_link_libraries(boundaryMesh ${ALL_LIBS})

add_executable(cube cube.cpp)
target_link_libraries(cube ${ALL_LIBS})

add_executable(damBreak damBreak.cpp)
target_link_libraries(damBreak ${ALL_LIBS})

add_executable(simpleBreakingDam_Ferstl2016 simpleBreakingDam_Ferstl2016.cpp)
target_link_libraries(simpleBreakingDam_Ferstl2016 ${ALL_LIBS})

add_executable(doubleDamBreak doubleDamBreak.cpp)
target_link_libraries(doubleDamBreak ${ALL_LIBS})

add_executable(sphere sphere.cpp)
target_link_libraries(sphere ${ALL_LIBS})

add_executable(monkey monkey.cpp)
target_link_libraries(monkey ${ALL_LIBS})

add_executable(takanaga takanaga.cpp)
target_link_libraries(takanaga ${ALL_LIBS})

add_executable(zeroGravity zeroGravity.cpp)
target_link_libraries(zeroGravity ${ALL_LIBS})

add_executable(twoZeroGravity twoZeroGravity.cpp)
target_link_libraries(twoZeroGravity ${ALL_LIBS})

add_executable(blenderExporter blenderExporter.cpp)
target_link_libraries(blenderExporter ${ALL_LIBS})

add_executable(hokusaiExporter hokusaiExporter.cpp)
target_link_libraries(hokusaiExporter ${ALL_LIBS})

add_executable(hokusaiImporter hokusaiImporter.cpp)
target_link_libraries(hokusaiImporter ${ALL_LIBS})

add_executable(debugExporter debugExporter.cpp)
target_link_libraries(debugExporter ${ALL_LIBS})

