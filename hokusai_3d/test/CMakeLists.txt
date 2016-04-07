cmake_minimum_required(VERSION 2.6)
project(test)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

#if(HOKUSAI_OPENMP)
#    find_package(OpenMP)
#    if(OPENMP_FOUND)
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    endif()
#endif()
#
#add_subdirectory(../hokusai hokusai/build)
#include_directories(../hokusai/include)
#
#find_package(Boost COMPONENTS timer system REQUIRED)
#include_directories( ${BOOST_INCLUDE_DIRS} )
#
#set(ALJABR_INCLUDE_DIRS "" CACHE PATH "path to alajabr headers")
#include_directories(${ALJABR_INCLUDE_DIRS})

add_subdirectory(gtest-1.7.0)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

set(
    ALL_LIBS
    #${Boost_LIBRARIES}
    #HOKUSAI 
    gtest
    gtest_main
    )

add_executable(vecTest vecTest.cpp)
target_link_libraries(vecTest ${ALL_LIBS})

add_test(VECTEST vecTest)
