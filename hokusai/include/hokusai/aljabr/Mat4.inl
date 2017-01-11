#ifndef ALJABR_MAT4_INL
#define ALJABR_MAT4_INL

#include "Mat4.hpp"
#include <stdexcept>

namespace aljabr
{

template<class T>
Mat4<T>::Mat4()
{
    m_data.fill(T(0));
}

template<class T>
Mat4<T>::Mat4(const Mat4<T>& m)
{
    m_data = m.m_data;
}

template<class T>
Mat4<T>::Mat4(const Mat3<T>& m, const unsigned int rowOffset, const unsigned int colOffset)
{
    try
    {
        identity();
        if(rowOffset>1 || colOffset>1) throw std::out_of_range("Bad offsets");
        for(size_t i=0; i<m.rows(); ++i)
            for(size_t j=0; j<m.cols(); ++j)
                m_data[index(rowOffset+i,colOffset+j)] = m(i,j);

    }
    catch(std::exception& e)
    {
        std::cout << "Standard exception: " << e.what() << std::endl;
    }
}

template<class T>
Mat4<T>::~Mat4()
{

}

template<class T>
T& Mat4<T>::operator()(unsigned int rowIndex, unsigned int colIndex)
{
    return m_data[index(rowIndex,colIndex)];
}

template<class T>
const T& Mat4<T>::operator()(unsigned int rowIndex, unsigned int colIndex) const
{
    return m_data[index(rowIndex,colIndex)];
}

template<class T>
T& Mat4<T>::operator[](unsigned int index)
{
    return m_data[index];
}

template<class T>
const T& Mat4<T>::operator[](unsigned int index) const
{
    return m_data[index];
}

template<class T>
Mat4<T>& Mat4<T>::operator=(const Mat4<T>& rhs)
{
    m_data = rhs.m_data;
    return *this;
}

template<class T>
Mat4<T>& Mat4<T>::operator*=(const Mat4<T>& rhs)
{
    Mat4<T> copy(*this);
    m_data.fill(0.0);
    size_t N = rows();
    for(size_t i=0; i<N; ++i)
        for(size_t j=0; j<N; ++j)
            for(size_t k=0; k<N; ++k)
                m_data[index(i,j)] += copy[index(i,k)]*rhs[index(k,j)];
    return *this;
}

template<class T>
size_t Mat4<T>::size() const
{
    return m_size;
}

template<class T>
size_t Mat4<T>::rows() const
{
    return m_rowsNumber;
}

template<class T>
size_t Mat4<T>::cols() const
{
    return m_colsNumber;
}

template<class T>
void Mat4<T>::fill(const T& value)
{
    m_data.fill(value);
}

template<class T>
unsigned int Mat4<T>::index(unsigned int rowIndex, unsigned int colIndex) const
{
    return rowIndex*4+colIndex;
}

template<class T>
void Mat4<T>::identity()
{
    for(size_t i=0; i<rows(); ++i)
        for(size_t j=0; j<cols(); ++j)
            if(i==j)
                m_data[index(i,j)] = T(1);
            else
                m_data[index(i,j)] = T(0);
}

template<class T>
void Mat4<T>::translation(const Vec3<T> &v)
{
    identity();
    for(size_t i=0; i<v.size(); ++i)
        m_data[index(i, cols()-1)] = v[i];
}

template<class T>
void Mat4<T>::rotation(const Vec3<T>& axis, const T& angle)
{
    Mat3<T> r3;
    r3.rotation(axis, angle);
    (*this) = Mat4<T>(r3);
}

template<class T>
Mat4<T> operator*(Mat4<T> lhs, const Mat4<T>& rhs)
{
    lhs *= rhs;
    return lhs;
}

template<class T>
Vec4<T> operator*(const Mat4<T>& lhs, const Vec4<T>& rhs)
{
    Vec4<T> r;
    r.fill(0.0);
    size_t N = rhs.size();
    for(size_t i=0; i<N; ++i)
        for(size_t j=0; j<N; ++j)
                r[i] += lhs(i,j)*rhs[j];
    return r;
}

template<class T>
Vec3<T> operator*(const Mat4<T>& lhs, const Vec3<T>& rhs)
{
    Vec4<T> r, rhs2(rhs[0], rhs[1], rhs[2], T(1.0) );
    r = lhs*rhs2;
    return Vec3<T>(r[0], r[1], r[2]);
}

}//namespace aljabr

#endif //ALJABR_MAT4_INL
