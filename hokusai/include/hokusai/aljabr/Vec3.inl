#ifndef ALJABR_VEC3_INL
#define ALJABR_VEC3_INL

#include "Vec3.hpp"

namespace aljabr
{

template<class T>
Vec3<T>::Vec3()
{
    m_data.fill(0);
}

template<class T>
Vec3<T>::Vec3(const T& x)
{
    m_data.fill(x);
}

template<class T>
Vec3<T>::Vec3(const Vec3& v)
{
    m_data = v.data();
}

template<class T>
Vec3<T>::Vec3(const T &x, const T &y, const T &z)
{
    m_data[0] = x;
    m_data[1] = y;
    m_data[2] = z;
}

template<class T>
Vec3<T>::~Vec3()
{}

template<class T>
size_t Vec3<T>::size() const
{
    return m_size;
}

template<class T>
void Vec3<T>::fill(const T& x)
{
    m_data.fill(x);
}

template<class T>
const std::array<T, 3>& Vec3<T>::data() const
{
    return m_data;
}

template<class T>
std::array<T, 3>& Vec3<T>::data()
{
    return m_data;
}

template<class T>
T Vec3<T>::lengthSquared() const
{
    return (m_data[0]*m_data[0]+m_data[1]*m_data[1]+m_data[2]*m_data[2]);
}

template<class T>
T Vec3<T>::length() const
{
    return std::sqrt(m_data[0]*m_data[0]+m_data[1]*m_data[1]+m_data[2]*m_data[2]);
}

template<class T>
void Vec3<T>::normalize()
{
    T l = std::sqrt(m_data[0]*m_data[0]+m_data[1]*m_data[1]+m_data[2]*m_data[2]); m_data[0]/=l; m_data[1]/=l; m_data[2]/=l;
}

template<class T>
Vec3<T>& Vec3<T>::operator= (const Vec3<T> & v)
{
    m_data = v.data();
    return *this;
}

template<class T>
Vec3<T>& Vec3<T>::operator*=(const T& factor)
{
    m_data[0]*=factor; m_data[1]*=factor; m_data[2]*=factor;
    return *this;
}

template<class T>
Vec3<T>& Vec3<T>::operator/=(const T& factor)
{
    m_data[0]/=factor; m_data[1]/=factor; m_data[2]/=factor;
    return *this;
}

template<class T>
Vec3<T>& Vec3<T>::operator+=(const Vec3<T> & v)
{
    m_data[0]+=v[0]; m_data[1]+=v[1]; m_data[2]+=v[2];
    return *this;
}

template<class T>
Vec3<T>& Vec3<T>::operator-=(const Vec3<T> & v)
{
    m_data[0]-=v[0]; m_data[1]-=v[1]; m_data[2]-=v[2];
    return *this;
}

template<class T>
T& Vec3<T>::operator[](int i)
{
    return m_data[i];
}

template<class T>
const T& Vec3<T>::operator[](int i) const
{
    return m_data[i];
}

template<class T>
T Vec3<T>::dotProduct ( const Vec3<T>& v1, const Vec3<T>& v2)
{
    return T(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] );
}

template<class T>
Vec3< T> Vec3<T>::crossProduct(const Vec3<T>& v1, const Vec3<T>& v2)
{
    return Vec3<T>(v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]);
}

template<class T>
Vec3<T> Vec3<T>::cwiseProd(const Vec3<T> &v) const
{
    return Vec3<T>( m_data[0]*v[0], m_data[1]*v[1], m_data[2]*v[2]);
}

template<class T>
Vec3<T> Vec3<T>::max(const T& value )
{
    return Vec3<T>( std::max( m_data[0], value ), std::max( m_data[1], value ), std::max( m_data[2], value ));
}

template<class T>
Vec3<T> Vec3<T>:: min(const T& value )
{
    return Vec3<T>( std::min( m_data[0], value ), std::min( m_data[1], value ), std::min( m_data[2], value ));
}

template<class T>
Vec3<T> Vec3<T>::abs( const Vec3<T>& v)
{
    return Vec3<T>(std::abs(v[0]), std::abs(v[1]), std::abs(v[2]));
}

template<class T>
Vec3<T> Vec3<T>::sign( const Vec3<T>& v)
{
    return Vec3<T>(signum<T>(v[0]), signum<T>(v[1]), signum<T>(v[2]));
}

template<class T>
T Vec3<T>::max(const Vec3<T>& v)
{
    return (*std::max_element(v.data().begin(), v.data().end()));
}

template<class T>
T Vec3<T>::min(const Vec3<T>& v)
{
    return (*std::min_element(v.data().begin(), v.data().end()));
}

//Non-member functions

template<class T>
const Vec3<T> operator* ( const T& factor, const Vec3<T> & vector )
{
    return Vec3<T>( factor*vector[0], factor*vector[1], factor*vector[2] );
}

template<class T>
const Vec3<T> operator* (const Vec3<T> & vector,  const T& factor )
{
    return Vec3<T>( factor*vector[0], factor*vector[1], factor*vector[2] );
}

template<class T>
const Vec3<T> operator- ( const Vec3<T> & vector )
{
    return Vec3< T>( -vector[0], -vector[1],-vector[2]);
}

template<class T>
const  Vec3<T> operator- ( const  Vec3<T> & v1, const  Vec3<T> & v2 )
{
    return Vec3<T>(v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]);
}

template<class T>
const  Vec3<T> operator+ ( const  Vec3<T> & v1, const  Vec3<T> & v2 )
{
    return Vec3<T>(v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]);
}


template<class T>
const Vec3<T> operator/ ( const Vec3<T> & vector, const T& factor)
{
    return Vec3<T>( vector[0]/factor, vector[1]/factor, vector[2]/factor );
}

template<class T>
bool operator== ( const Vec3<T>& v1, const Vec3<T>& v2 )
{
    return ( ( v1[0]==v2[0] ) && ( v1[1]==v2[1]) && ( v1[2]==v2[2]) );
}

template<class T>
bool operator!= ( const Vec3<T>& v1, const Vec3<T>& v2 )
{
    return ( ( v1[0]!=v2[0] ) && ( v1[1]!=v2[1]) && ( v1[2]!=v2[2]) );
}

template<class T>
std::istream& operator >> ( std::istream& in, Vec3<T>& v )
{
    for( int i=0; i<3; ++i ) in>>v[i];
    return in;
}

template<class T>
std::ostream& operator << ( std::ostream& out, const Vec3<T>& v )
{
    out << v[0] << " " << v[1] << " " << v[2];
    return out;
}

}//namespace aljabr

#endif //ALJABR_VEC3_INL
