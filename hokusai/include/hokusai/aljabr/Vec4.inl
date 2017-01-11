#ifndef ALJABR_VEC4_INL
#define ALJABR_VEC4_INL

#include "Vec4.hpp"

namespace aljabr
{

template<class T>
Vec4<T>::Vec4()
{
    m_data.fill(T(0));
}

template<class T>
Vec4<T>::Vec4(const T& x, const T& y, const T& z, const T& w)
{
    m_data[0] = x;
    m_data[1] = y;
    m_data[2] = z;
    m_data[3] = w;
}

template<class T>
Vec4<T>::~Vec4()
{

}

template<class T>
T& Vec4<T>::operator[](unsigned int index)
{
    return m_data[index];
}

template<class T>
const T& Vec4<T>::operator[](unsigned int index) const
{
    return m_data[index];
}

template<class T>
size_t Vec4<T>::size() const
{
    return m_size;
}

template<class T>
void Vec4<T>::fill(const T& value)
{
    m_data.fill(value);
}

}//namespace aljabr

#endif //ALJABR_VEC4_INL
