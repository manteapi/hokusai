#ifndef ALJABR_VEC2_HPP
#define ALJABR_VEC2_HPP

#include <cmath>
#include <iostream>
#include <array>
#include <algorithm>
#include "Utils.hpp"

namespace aljabr
{

template<class T>
class Vec2
{
private:
    static const size_t m_size = 2;
    std::array<T, m_size> m_data;

public:
    Vec2();
    Vec2(const T& x);
    Vec2(const Vec2& v);
    Vec2(const T& x, const T& y);
    ~Vec2();

    size_t size() const;
    void fill(const T& x);

    const std::array<T, 2>& data() const;
    std::array<T, 2>& data();

    T lengthSquared() const;
    T length() const;
    void normalize();

    Vec2<T> cwiseProd(const Vec2<T> &v) const;
    Vec2<T> max(const T& value );
    Vec2<T> min(const T& value );

    Vec2<T>& operator= (const Vec2<T> & v);
    Vec2<T>& operator*=(const T& factor);
    Vec2<T>& operator/=(const T& factor);
    Vec2<T>& operator+=(const Vec2<T> & v);
    Vec2<T>& operator-=(const Vec2<T> & v);
    bool operator< (const Vec2<T>& v ) const;
    T& operator[](int i);
    const T& operator[](int i) const;

    static T dotProduct ( const Vec2<T>& v1, const Vec2<T>& v2);
    static T crossProduct(const Vec2<T>& v1, const Vec2<T>& v2);
    static T max(const Vec2<T>& v);
    static T min(const Vec2<T>& v);
    static Vec2<T> abs(const Vec2<T>& v);
    static Vec2<T> sign(const Vec2<T>& v);
};

template<class T>
const Vec2<T> operator* (const T& factor, const Vec2<T> & vector );

template<class T>
const Vec2<T> operator* (const Vec2<T> & vector, const T& factor );

template<class T>
const Vec2<T> operator- ( const Vec2<T> & vector );

template<class T>
const  Vec2<T> operator- ( const  Vec2<T> & v1, const  Vec2<T> & v2 );

template<class T>
const Vec2<T> operator/ ( const Vec2<T> & vector, const T& factor );

template<class T>
bool operator== ( const Vec2<T>& v1, const Vec2<T>& v2 );

template<class T>
bool operator!= ( const Vec2<T>& v1, const Vec2<T>& v2 );

template<class T>
std::istream& operator >> ( std::istream& in, Vec2<T>& v );

template<class T>
std::ostream& operator << ( std::ostream& out, const Vec2<T>& v );


}//aljabr

#endif // ALJABR_VEC2_HPP
