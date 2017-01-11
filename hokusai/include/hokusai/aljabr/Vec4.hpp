#ifndef ALJABR_VEC4_HPP
#define ALJABR_VEC4_HPP

#include <cmath>
#include <iostream>
#include <array>
#include <algorithm>

namespace aljabr
{

template<class T>
class Vec4
{
private:
    static const size_t m_size = 4;
    std::array< T, m_size > m_data;

public:
    Vec4();
    ~Vec4();
    Vec4(const T &x, const T &y, const T &z, const T &w);

    size_t size() const;

    void fill(const T& value);

    T& operator[](unsigned int index);
    const T& operator[](unsigned int index) const;
};

}//aljabr

#endif // ALJABR_VEC4_HPP
