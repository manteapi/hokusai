cmake_minimum_required(VERSION 2.6)
project(fluid)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-std=c++0x -fopenmp" ${CMAKE_CXX_FLAGS})

add_subdirectory(hokusai hokusai/build)
include_directories(hokusai/include)

#find_package(OpenGL REQUIRED)
#include_directories( ${OPENGL_INCLUDE_DIRS} )

find_package(Boost COMPONENTS timer system REQUIRED)
include_directories( ${BOOST_INCLUDE_DIRS} )

set(
    ALL_LIBS
    ${Boost_LIBRARIES}
    HOKUSAI 
    )

add_executable(damBreak examples/damBreak.cpp)
target_link_libraries(damBreak ${ALL_LIBS})

add_executable(cube examples/cube.cpp)
target_link_libraries(cube ${ALL_LIBS})

add_executable(takanaga examples/takanaga.cpp)
target_link_libraries(takanaga ${ALL_LIBS})

add_executable(zeroGravity examples/zeroGravity.cpp)
target_link_libraries(zeroGravity ${ALL_LIBS})

add_executable(monkey examples/monkey.cpp)
target_link_libraries(monkey ${ALL_LIBS})
