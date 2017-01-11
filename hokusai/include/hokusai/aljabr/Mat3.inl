#ifndef ALJABR_MAT3_INL
#define ALJABR_MAT3_INL

#include "Mat3.hpp"

namespace aljabr
{

template<class T>
Mat3<T>::Mat3()
{
    m_data.fill(T(0));
}

template<class T>
Mat3<T>::~Mat3()
{

}

template<class T>
Mat3<T>::Mat3(const Mat3& m)
{
    m_data = m.m_data;
}

template<class T>
T& Mat3<T>::operator()(unsigned int rowIndex, unsigned int colIndex)
{
    return m_data[index(rowIndex,colIndex)];
}

template<class T>
const T& Mat3<T>::operator()(unsigned int rowIndex, unsigned int colIndex) const
{
    return m_data[index(rowIndex,colIndex)];
}

template<class T>
T& Mat3<T>::operator[](unsigned int index)
{
    return m_data[index];
}

template<class T>
const T& Mat3<T>::operator[](unsigned int index) const
{
    return m_data[index];
}

template<class T>
Mat3<T>& Mat3<T>::operator=(const Mat3<T>& rhs)
{
    m_data = rhs.m_data;
    return *this;
}

template<class T>
Mat3<T>& Mat3<T>::operator*=(const T& s)
{
    for(size_t i=0; i<m_size; ++i) m_data[i] *= s;
    return *this;
}

template<class T>
Mat3<T>& Mat3<T>::operator+=(const Mat3<T>& rhs)
{
    for(size_t i=0; i<m_size; ++i) m_data[i] += rhs[i];
    return *this;
}

template<class T>
size_t Mat3<T>::size() const
{
    return m_size;
}

template<class T>
size_t Mat3<T>::rows() const
{
    return m_rowsNumber;
}

template<class T>
size_t Mat3<T>::cols() const
{
    return m_colsNumber;
}

template<class T>
unsigned int Mat3<T>::index(unsigned int rowIndex, unsigned int colIndex) const
{
    return rowIndex*3+colIndex;
}

template<class T>
void Mat3<T>::identity()
{
    for(size_t i=0; i<rows(); ++i)
        for(size_t j=0; j<cols(); ++j)
            if(i==j)
                m_data[index(i,j)] = T(1);
            else
                m_data[index(i,j)] = T(0);
}

template<class T>
void Mat3<T>::cross(const Vec3<T>& v)
{
    for(size_t i=0; i<rows(); ++i) m_data[index(i,i)] = T(0);
    m_data[index(0,1)] = -v[2];
    m_data[index(0,2)] = v[1];
    m_data[index(1,0)] = v[2];
    m_data[index(1,2)] = -v[0];
    m_data[index(2,0)] = -v[1];
    m_data[index(2,1)] = v[0];
}

template<class T>
void Mat3<T>::tensor(const Vec3<T>& v)
{
    for(size_t i=0; i<rows(); ++i)
        for(size_t j=0; j<cols(); ++j)
            m_data[index(i,j)] = v[i]*v[j];
}

template<class T>
void Mat3<T>::rotation(const Vec3<T>& axis, const T& theta)
{
    Mat3<T> i;
    i.identity();

    Mat3<T> c;
    c.cross(axis);

    Mat3<T> t;
    t.tensor(axis);

    Mat3<T> r;
    r += std::cos( (T)theta)*i;
    r += std::sin( (T)theta )*c;
    r += ( (T)(1.0) - std::cos( (T)(theta) ) )*t;

    (*this) = r;
}

template<class T>
Mat3<T> operator*(Mat3<T> lhs, const T& rhs)
{
    lhs*=rhs;
    return lhs;
}

template<class T>
Mat3<T> operator*(const T& lhs, Mat3<T> rhs)
{
    rhs*=lhs;
    return rhs;
}

template<class T>
Mat3<T> operator+(Mat3<T> lhs, const Mat3<T>& rhs)
{
    lhs+=rhs;
    return lhs;
}

}//namespace aljabr

#endif //ALJABR_MAT3_INL
