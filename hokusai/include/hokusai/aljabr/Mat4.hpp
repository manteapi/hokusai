#ifndef ALJABR_MAT4_HPP
#define ALJABR_MAT4_HPP

#include <cmath>
#include <iostream>
#include <array>
#include <algorithm>
#include "Vec3.inl"
#include "Vec4.inl"
#include "Mat3.inl"

namespace aljabr
{

///Ordering
///     | 0  1  2  3  |
///     | 4  5  6  7  |
/// M = | 8  9  10 11 |
///     | 12 13 14 15 |
///
/// index = rowIndex*4+columnIndex (Row major)

template<class T>
class Mat4
{
private:
    static const size_t m_rowsNumber = 4;
    static const size_t m_colsNumber = 4;
    static const size_t m_size = 16;
    std::array< T, m_size > m_data;

public:
    Mat4();
    Mat4(const Mat4<T>& m);
    Mat4(const Mat3<T>& m, const unsigned int rowOffset=0, const unsigned int colOffset=0);
    ~Mat4();

    size_t size() const;
    size_t rows() const;
    size_t cols() const;

    void fill(const T& value);
    unsigned int index(unsigned int rowIndex, unsigned int colIndex) const;

    T& operator[](unsigned int index);
    const T& operator[](unsigned int index) const;

    T& operator()(unsigned int rowIndex, unsigned int colIndex);
    const T& operator()(unsigned int rowIndex, unsigned int colIndex) const;

    Mat4<T>& operator=(const Mat4<T>& rhs);
    Mat4<T>& operator*=(const Mat4<T>& rhs);

    void identity();
    void translation(const Vec3<T>& v);
    void rotation(const Vec3<T>& axis, const T& angle);
};

template<class T>
Mat4<T> operator*(Mat4<T> lhs, const Mat4<T>& rhs);
template<class T>
Vec4<T> operator*(const Mat4<T>& lhs, const Vec4<T>& rhs);
template<class T>
Vec3<T> operator*(const Mat4<T>& lhs, const Vec3<T>& rhs);

}//aljabr

#endif // ALJABR_MAT4_HPP
