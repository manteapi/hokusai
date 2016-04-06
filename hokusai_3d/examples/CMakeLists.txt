cmake_minimum_required(VERSION 2.6)
project(hokusai_examples)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

add_subdirectory(../hokusai hokusai/build)
include_directories(../hokusai/include)

find_package(Boost COMPONENTS timer system REQUIRED)
include_directories( ${BOOST_INCLUDE_DIRS} )

set(ALJABR_INCLUDE_DIRS "" CACHE PATH "path to alajabr headers")
include_directories(${ALJABR_INCLUDE_DIRS})

set(MAGNET_INCLUDE_DIRS "" CACHE PATH "path to mosaic headers")
include_directories(${MAGNET_INCLUDE_DIRS})

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

add_executable(particleSource particleSource.cpp)
target_link_libraries(particleSource ${ALL_LIBS})

add_executable(boundaryMesh boundaryMesh.cpp)
target_link_libraries(boundaryMesh ${ALL_LIBS})

add_executable(cube cube.cpp)
target_link_libraries(cube ${ALL_LIBS})

add_executable(damBreak damBreak.cpp)
target_link_libraries(damBreak ${ALL_LIBS})

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
