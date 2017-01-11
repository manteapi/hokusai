#ifndef ALJABR_VEC3_HPP
#define ALJABR_VEC3_HPP

#include <cmath>
#include <iostream>
#include <array>
#include <algorithm>
#include "Utils.hpp"

namespace aljabr
{

template<class T>
class Vec3
{
private:
    static const size_t m_size = 3;
    std::array<T, m_size > m_data;

public:
    Vec3();
    Vec3(const T& x);
    Vec3(const Vec3& v);
    Vec3(const T& x, const T& y, const T& z);
    ~Vec3();

    size_t size() const;
    void fill(const T& x);

    const std::array<T, 3>& data() const;
    std::array<T, 3>& data();

    T lengthSquared() const;
    T length() const;
    void normalize();

    Vec3<T> cwiseProd(const Vec3<T> &v) const;
    Vec3<T> max(const T& value);
    Vec3<T> min(const T& value);

    Vec3<T>& operator= (const Vec3<T>& v);
    Vec3<T>& operator*=(const T& factor);
    Vec3<T>& operator/=(const T& factor);
    Vec3<T>& operator+=(const Vec3<T> & v);
    Vec3<T>& operator-=(const Vec3<T> & v);
    bool operator < (const Vec3<T>& v ) const;
    T& operator[](int i);
    const T& operator[](int i) const;

    static T dotProduct ( const Vec3<T>& v1, const Vec3<T>& v2);
    static Vec3<T> crossProduct(const Vec3<T>& v1, const Vec3<T>& v2);
    static T max(const Vec3<T>& v);
    static T min(const Vec3<T>& v);
    static Vec3<T> abs( const Vec3<T>& v);
    static Vec3<T> sign( const Vec3<T>& v);
};

template<class T>
const Vec3<T> operator* (const T& factor, const Vec3<T> & vector );
template<class T>
const Vec3<T> operator* (const Vec3<T> & vector, const T& factor );
template<class T>
const Vec3<T> operator- ( const Vec3<T> & vector );
template<class T>
const  Vec3<T> operator- ( const Vec3<T> & v1, const Vec3<T> & v2 );
template<class T>
const  Vec3<T> operator+ ( const Vec3<T> & v1, const Vec3<T> & v2 );
template<class T>
const Vec3<T> operator/ ( const Vec3<T> & vector, const T& factor );
template<class T>
bool operator== ( const Vec3<T>& v1, const Vec3<T>& v2 );
template<class T>
bool operator!= ( const Vec3<T>& v1, const Vec3<T>& v2 );
template<class T>
std::istream& operator >> ( std::istream& in, Vec3<T>& v );
template<class T>
std::ostream& operator << ( std::ostream& out, const Vec3<T>& v );

}//aljabr

#endif // ALJABR_VEC3_HPP
