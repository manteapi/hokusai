cmake_minimum_required(VERSION 2.6)
project(test)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_subdirectory(../hokusai hokusai/build)
include_directories(../hokusai/include)

add_subdirectory(gtest-1.7.0)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

set(
    ALL_LIBS
    ${Boost_LIBRARIES}
    HOKUSAI 
    gtest
    gtest_main
    )

add_executable(test1 test1.cpp)
target_link_libraries(test1 ${ALL_LIBS})

add_executable(test2 test2.cpp)
target_link_libraries(test2 ${ALL_LIBS})

add_test(TEST1 test1)
add_test(TEST2 test2)
