#ifndef ALJABR_MAT3_HPP
#define ALJABR_MAT3_HPP

#include <cmath>
#include <iostream>
#include <array>
#include <algorithm>
#include "Vec3.hpp"

namespace aljabr
{

///Ordering
///     | 0  1  2 |
/// M = | 3  4  5 |
///     | 6  7  8 |
/// index = rowIndex*3+columnIndex (Row major)

template<class T>
class Mat3
{
private:
    static const size_t m_rowsNumber = 3;
    static const size_t m_colsNumber = 3;
    static const size_t m_size = 9;
    std::array< T, m_size > m_data;

public:
    Mat3();
    Mat3(const Mat3& m);
    ~Mat3();

    size_t size() const;
    size_t rows() const;
    size_t cols() const;

    unsigned int index(unsigned int rowIndex, unsigned int colIndex) const;

    T& operator[](unsigned int index);
    const T& operator[](unsigned int index) const;

    T& operator()(unsigned int rowIndex, unsigned int colIndex);
    const T& operator()(unsigned int rowIndex, unsigned int colIndex) const;

    Mat3<T>& operator=(const Mat3<T>& rhs);
    Mat3<T>& operator*=(const T& s);
    Mat3<T>& operator+=(const Mat3<T>& rhs);

    void identity();
    void cross(const Vec3<T>& v);
    void tensor(const Vec3<T>& v);
    void rotation(const Vec3<T>& axis, const T &theta);

};

template<class T>
Mat3<T> operator*(Mat3<T> lhs, const T& rhs);
template<class T>
Mat3<T> operator*(const T& lhs, Mat3<T> rhs);

template<class T>
Mat3<T> operator+(Mat3<T> lhs, const Mat3<T>& rhs);

}//aljabr

#endif // ALJABR_MAT3_HPP
