cmake_minimum_required(VERSION 2.6)
project(online_examples)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

set(MAGNET_INCLUDE_DIRS "" CACHE PATH "path to magnet headers")
include_directories(${MAGNET_INCLUDE_DIRS})

set(EIGEN_INCLUDE_DIRS "" CACHE PATH "path to eigen headers")
include_directories(${EIGEN_INCLUDE_DIRS})

set(HOKUSAI_INCLUDE_DIRS "" CACHE PATH "path to hokusai headers")
include_directories(${HOKUSAI_INCLUDE_DIRS})
set(HOKUSAI_LIBRARIES "" CACHE FILEPATH "path to hokusailibrary")

set(GLM_INCLUDE_DIRS "" CACHE PATH "path to glm headers")
include_directories(${GLM_INCLUDE_DIRS})

set(GLEW_INCLUDE_DIRS "" CACHE PATH "path to glew headers")
set(GLEW_LIBRARIES "" CACHE FILEPATH "path to glew library")
include_directories(${GLEW_INCLUDE_DIRS})

set(SFML_INCLUDE_DIRS "" CACHE PATH "path to sfml headers")
set(SFML_SYSTEM_LIBRARIES "" CACHE FILEPATH "path to sfml system library")
set(SFML_WINDOW_LIBRARIES "" CACHE FILEPATH "path to sfml windo library")
set(SFML_GRAPHICS_LIBRARIES "" CACHE FILEPATH "path to sfml graphics library")
include_directories(${SFML_INCLUDE_DIRS})

set(SOURCES
    main.cpp
)

set(HEADERS
    shader.hpp
    shader.cpp
)

set(SHADERS
    particle_vertex.glsl
    particle_fragment.glsl
    grid_vertex.glsl
    grid_fragment.glsl
)

set(EXECUTABLE_NAME main)
add_executable(${EXECUTABLE_NAME} ${SOURCES} ${HEADERS} ${SHADERS})

find_package(OpenGL REQUIRED)
if (OPENGL_FOUND)
    target_link_libraries(${EXECUTABLE_NAME} ${OPENGL_LIBRARIES})
    target_link_libraries(${EXECUTABLE_NAME} m)  # if you use maths.h
endif()

if(HOKUSAI_OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

target_link_libraries(${EXECUTABLE_NAME} ${SFML_SYSTEM_LIBRARIES})
target_link_libraries(${EXECUTABLE_NAME} ${SFML_WINDOW_LIBRARIES})
target_link_libraries(${EXECUTABLE_NAME} ${SFML_GRAPHICS_LIBRARIES})
target_link_libraries(${EXECUTABLE_NAME} ${GLEW_LIBRARIES})
target_link_libraries(${EXECUTABLE_NAME} ${HOKUSAI_LIBRARIES})

