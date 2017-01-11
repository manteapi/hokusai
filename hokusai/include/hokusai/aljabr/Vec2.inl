#ifndef ALJABR_VEC2_INL
#define ALJABR_VEC2_INL

#include "Vec2.hpp"

namespace aljabr
{

template<class T>
Vec2<T>::Vec2()
{
    m_data.fill(0);
}

template<class T>
Vec2<T>::Vec2(const T& x)
{
    m_data.fill(x);
}

template<class T>
Vec2<T>::Vec2(const Vec2& v)
{
    m_data = v.data();
}

template<class T>
Vec2<T>::Vec2(const T &x, const T &y)
{
    m_data[0] = x; m_data[1] = y;
}

template<class T>
Vec2<T>::~Vec2()
{}

template<class T>
size_t Vec2<T>::size() const
{
    return m_size;
}

template<class T>
void Vec2<T>::fill(const T& x)
{
    m_data.fill(x);
}

template<class T>
const std::array<T, 2>& Vec2<T>::data() const
{
    return m_data;
}

template<class T>
std::array<T, 2>& Vec2<T>::data()
{
    return m_data;
}

template<class T>
T Vec2<T>::lengthSquared() const
{
    return (m_data[0]*m_data[0]+m_data[1]*m_data[1]);
}

template<class T>
T Vec2<T>::length() const
{
    return std::sqrt(m_data[0]*m_data[0]+m_data[1]*m_data[1]);
}

template<class T>
void Vec2<T>::normalize()
{
    T l = std::sqrt(m_data[0]*m_data[0]+m_data[1]*m_data[1]);
    (*this)/=l;
}

template<class T>
Vec2<T>& Vec2<T>::operator= (const Vec2<T> & v)
{
    m_data = v.data();
    return *this;
}

template<class T>
Vec2<T>& Vec2<T>::operator*=(const T& factor)
{
    m_data[0]*=factor; m_data[1]*=factor;
    return *this;
}

template<class T>
Vec2<T>& Vec2<T>::operator/=(const T& factor)
{
    m_data[0]/=factor; m_data[1]/=factor;
    return *this;
}

template<class T>
Vec2<T>& Vec2<T>::operator+=(const Vec2<T>& v)
{
    m_data[0]+=v[0]; m_data[1]+=v[1];
    return *this;
}

template<class T>
Vec2<T>& Vec2<T>::operator-=(const Vec2<T>& v)
{
    m_data[0]-=v[0]; m_data[1]-=v[1];
    return *this;
}

template<class T>
bool Vec2<T>::operator< (const Vec2<T>& v ) const
{
    for( int i=0; i<2; i++ )
    {
        if( m_data[i]<v[i] )
            return true;
        else if( m_data[i]>v[i] )
            return false;
    }
    return false;
}

template<class T>
T& Vec2<T>::operator[](int i)
{
    return m_data[i];
}

template<class T>
const T& Vec2<T>::operator[](int i) const
{
    return m_data[i];
}

template<class T>
T Vec2<T>::dotProduct ( const Vec2<T>& v1, const Vec2<T>& v2)
{
    return( v1[0]*v2[0]+v1[1]*v2[1]);
}

template<class T>
T Vec2<T>::crossProduct(const Vec2<T>& v1, const Vec2<T>& v2)
{
    return (v1[0]*v2[1] - v1[1]*v2[0]);
}

template<class T>
inline Vec2<T> Vec2<T>::cwiseProd(const Vec2<T> &v) const
{
    return Vec2<T>( m_data[0]*v[0], m_data[1]*v[1]);
}

template<class T>
inline Vec2<T> Vec2<T>::max(const T& value )
{
    return Vec2<T>( std::max( m_data[0], value ), std::max( m_data[1], value ));
}

template<class T>
inline Vec2<T> Vec2<T>::min(const T& value )
{
    return Vec2<T>( std::min( m_data[0], value ), std::min( m_data[1], value ));
}

template<class T>
T Vec2<T>::max(const Vec2<T>& v)
{
    return (*std::max_element(v.data().begin(), v.data().end()));
}

template<class T>
T Vec2<T>::min(const Vec2<T>& v)
{
    return (*std::min_element(v.data().begin(), v.data().end()));
}

template<class T>
Vec2<T> Vec2<T>::abs( const Vec2<T>& v)
{
    return Vec2<T>(std::abs(v[0]), std::abs(v[1]));
}

template<class T>
Vec2<T> Vec2<T>::sign( const Vec2<T>& v)
{
    return Vec2<T>(signum<T>(v[0]), signum<T>(v[1]));
}

//Non-member functions

template<class T>
const Vec2<T> operator* (const T& factor, const Vec2<T> & vector )
{
    return Vec2<T>( factor*vector[0], factor*vector[1] );
}

template<class T>
const Vec2<T> operator* (const Vec2<T> & vector,  const T& factor)
{
    return Vec2<T>( factor*vector[0], factor*vector[1] );
}

template<class T>
const Vec2<T> operator- ( const Vec2<T> & vector )
{
    return Vec2< T>( -vector[0], -vector[1] );
}

template<class T>
const  Vec2<T> operator- ( const  Vec2<T> & v1, const  Vec2<T> & v2 )
{
    return Vec2<T>(v1[0]-v2[0], v1[1]-v2[1]);
}

template<class T>
const Vec2<T> operator/ ( const Vec2<T> & vector, const T& factor )
{
    return Vec2<T>( vector[0]/factor, vector[1]/factor );
}

template<class T>
bool operator== ( const Vec2<T>& v1, const Vec2<T>& v2 )
{
    return ( ( v1[0]==v2[0] ) && ( v1[1]==v2[1]) );
}

template<class T>
bool operator!= ( const Vec2<T>& v1, const Vec2<T>& v2 )
{
    return ( ( v1[0]!=v2[0] ) && ( v1[1]!=v2[1]) );
}

template<class T>
std::istream& operator >> ( std::istream& in, Vec2<T>& v )
{
    in>>v[0];
    in>>v[1];
    return in;
}

template<class T>
std::ostream& operator << ( std::ostream& out, const Vec2<T>& v )
{
    out << v[0] << " " << v[1];
    return out;
}

}//namespace aljabr

#endif //ALJABR_VEC2_INL
