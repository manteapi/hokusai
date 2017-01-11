cmake_minimum_required(VERSION 2.8)

project(aljabr)

set(CMAKE_CXX_FLAGS "-std=c++0x")

set(ALJABR_HEADER
    Vec2.hpp
    Vec2.inl
    Vec3.hpp
    Vec3.inl
    Vec4.hpp
    Vec4.inl
    Vector.hpp
    Vector.inl
    Mat2.hpp
    Mat2.inl
    Mat3.hpp
    Mat3.inl
    Mat4.hpp
    Mat4.inl
    Matrix.hpp
    Matrix.inl
    Utils.hpp
    Utils.inl
    AljabrCore
)

add_library(ALJABR_LIBRARY STATIC ${ALJABR_HEADER})
set_target_properties(ALJABR_LIBRARY PROPERTIES LINKER_LANGUAGE CXX)
