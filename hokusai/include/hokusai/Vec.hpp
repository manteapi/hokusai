#ifndef VEC_H
#define VEC_H

#include <cmath>
#include <iostream>
#include <array>
#include <algorithm>

using namespace std;

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//						Vec3< Real >
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

template<class Real>
int signum(Real val)
{
    return (Real(0)<val)-(val<Real(0));
}

template<class Real>
class Vec3
{
private:
    array< Real, 3 > vector;

public:
    inline Vec3(){ vector.fill(0);}
    inline Vec3(Real x){vector.fill(x);}
    inline Vec3(const Vec3& v){ vector = v.getConstArray();  }
    inline Vec3( Real x, Real y, Real z){ vector[0] = x; vector[1] = y; vector[2] = z; }
    ~Vec3(){}

    inline void info() const{ std::cout << "Vec : " << vector[0] << ", " << vector[1] << ", " << vector[2] << std::endl; }

    inline const array<Real, 3>& getConstArray() const { return vector; }
    inline array<Real, 3>& getArrayValue(){ return vector; }
    inline Real lengthSquared() const { return (vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);}
    inline Real length() const { return sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);}
    inline void normalize(){ Real l = sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]); vector[0]/=l; vector[1]/=l; vector[2]/=l;}
    inline void setAllValue(Real x){ vector[0] = vector[1] = vector[2] = x; }

    inline Vec3< Real >& operator= (const Vec3<Real> & v){ vector = v.getConstArray(); return *this;}

    inline Vec3< Real > operator-(const Vec3<Real> & v) const{ return Vec3<Real>( vector[0]-v[0], vector[1]-v[1], vector[2]-v[2] ); }
    inline Vec3< Real > operator+(const Vec3<Real> & v) const{ return Vec3<Real>( vector[0]+v[0], vector[1]+v[1], vector[2]+v[2] ); }

    inline Vec3< Real > operator*=(Real factor) { vector[0]*=factor; vector[1]*=factor; vector[2]*=factor; return *this; }
    inline Vec3< Real > operator/=(Real factor) { vector[0]/=factor; vector[1]/=factor; vector[2]/=factor; return *this; }

    inline Vec3< Real > operator+=(const Vec3<Real> & v) { vector[0]+=v[0]; vector[1]+=v[1]; vector[2]+=v[2]; return *this; }
    inline Vec3< Real > operator-=(const Vec3<Real> & v) { vector[0]-=v[0]; vector[1]-=v[1]; vector[2]-=v[2]; return *this; }

    inline Real& operator[](int i){ return vector[i]; }
    inline const Real& operator[](int i) const{ return vector[i]; }

    inline static Real dotProduct ( const Vec3<Real>& v1, const Vec3<Real>& v2){return Real(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] );}

    inline Vec3< Real > cwiseProd(const Vec3<Real> &v) const { return Vec3<Real>( vector[0]*v[0], vector[1]*v[1], vector[2]*v[2]);}

    inline Vec3< Real > max( Real value ){ return Vec3<Real>( std::max( vector[0], value ), std::max( vector[1], value ), std::max( vector[2], value ));}

    inline Vec3< Real > min( Real value ){ return Vec3<Real>( std::min( vector[0], value ), std::min( vector[1], value ), std::min( vector[2], value ));}

    inline static Vec3< Real > abs( const Vec3<Real>& v){return Vec3<Real>(std::abs(v[0]), std::abs(v[1]), std::abs(v[2]));}
    inline static Vec3< Real > sign( const Vec3<Real>& v){return Vec3<Real>(signum<Real>(v[0]), signum<Real>(v[1]), signum<Real>(v[2]));}

    inline static Real max(Vec3<Real> v){ return (*std::max_element(v.getConstArray().begin(), v.getConstArray().end())); }
    inline static Real min(Vec3<Real> v){ return (*std::min_element(v.getConstArray().begin(), v.getConstArray().end())); }
};

template<class Real>
inline const Vec3<Real>	operator* ( Real factor, const Vec3< Real > & vector ){return Vec3<Real>( factor*vector[0], factor*vector[1], factor*vector[2] );}
template<class Real>
inline const Vec3<Real>	operator* (const Vec3< Real > & vector,  Real factor ){return Vec3<Real>( factor*vector[0], factor*vector[1], factor*vector[2] );}
template<class Real>
inline const Vec3<Real>	operator- ( const Vec3<Real> & vector ){ return Vec3< Real>( -vector[0], -vector[1],-vector[2]); }
template<class Real>
inline const  Vec3<Real> operator- ( const  Vec3<Real> & v1, const  Vec3<Real> & v2 ){ return Vec3< Real >(v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]); }
template<class Real>
inline const Vec3<Real>	operator/ ( const Vec3<Real> & vector, Real divisor ){ return Vec3<Real>( vector[0]/divisor, vector[1]/divisor, vector[2]/divisor ); }
template<class Real> inline bool operator== ( const Vec3< Real >& v1, const Vec3< Real >& v2 ){ return ( ( v1[0]==v2[0] ) && ( v1[1]==v2[1]) && ( v1[2]==v2[2]) ); }
template<class Real> inline bool operator!= ( const Vec3< Real >& v1, const Vec3< Real >& v2 ){ return ( ( v1[0]!=v2[0] ) && ( v1[1]!=v2[1]) && ( v1[2]!=v2[2]) ); }

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//						Vec2< Real >
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

template<class Real>
class Vec2
{
private:
    array<Real, 2> vector;

public:
    inline Vec2(){vector.fill(0);}
    inline Vec2(Real x){vector.fill(x);}
    inline Vec2(const Vec2& v){ vector = v.getConstArray();  }
    inline Vec2( Real x, Real y){ vector[0] = x; vector[1] = y; }
    ~Vec2(){}

    inline void info(){ std::cout << "Vec : " << vector[0] << ", " << vector[1] << std::endl; }

    inline const array<Real, 2>& getConstArray() const { return vector; }
    inline array<Real, 2>& getArrayValue(){ return vector; }
    inline Real lengthSquared() const { return (vector[0]*vector[0]+vector[1]*vector[1]);}
    inline Real length() const { return sqrt(vector[0]*vector[0]+vector[1]*vector[1]);}
    inline void normalize(){ Real l = sqrt(vector[0]*vector[0]+vector[1]*vector[1]); vector[0]/=l; vector[1]/=l;}
    inline void setAllValue(Real x){ vector[0] = vector[1] = x; }

    inline Vec2< Real >& operator= (const Vec2<Real> & v){ vector = v.getConstArray(); return *this;}

    inline Vec2< Real > operator-(const Vec2<Real> & v) const{ return Vec2<Real>( vector[0]-v[0], vector[1]-v[1] ); }
    inline Vec2< Real > operator+(const Vec2<Real> & v) const{ return Vec2<Real>( vector[0]+v[0], vector[1]+v[1] ); }

    inline Vec2< Real > operator*=(Real factor) { vector[0]*=factor; vector[1]*=factor; return *this; }
    inline Vec2< Real > operator/=(Real factor) { vector[0]/=factor; vector[1]/=factor; return *this; }

    inline Vec2< Real > operator+=(const Vec2<Real> & v) { vector[0]+=v[0]; vector[1]+=v[1]; return *this; }
    inline Vec2< Real > operator-=(const Vec2<Real> & v) { vector[0]-=v[0]; vector[1]-=v[1]; return *this; }

    inline Real& operator[](int i){ return vector[i]; }
    inline const Real& operator[](int i) const{ return vector[i]; }

    inline static Real dotProduct ( const Vec2<Real>& v1, const Vec2<Real>& v2){return( v1[0]*v2[0]+v1[1]*v2[1]);}
    inline Vec2< Real > cwiseProd(const Vec2<Real> &v) const { return Vec2<Real>( vector[0]*v[0], vector[1]*v[1]);}
    inline Vec2< Real > maxVector( Real value ){ return Vec2<Real>( std::max( vector[0], value ), std::max( vector[1], value ));}
    inline Vec2< Real > minVector( Real value ){ return Vec2<Real>( std::min( vector[0], value ), std::min( vector[1], value ));}

    inline static Real max(Vec2<Real> v){ return (*std::max_element(v.getConstArray().begin(), v.getConstArray().end())); }
    inline static Real min(Vec2<Real> v){ return (*std::min_element(v.getConstArray().begin(), v.getConstArray().end())); }
    inline static Vec2< Real > abs( const Vec2<Real>& v){return Vec2<Real>(std::abs(v[0]), std::abs(v[1]));}
    inline static Vec2< Real > sign( const Vec2<Real>& v){return Vec2<Real>(signum<Real>(v[0]), signum<Real>(v[1]));}
};

template<class Real>
inline const Vec2<Real>	operator* ( Real factor, const Vec2< Real > & vector ){return Vec2<Real>( factor*vector[0], factor*vector[1] );}
template<class Real>
inline const Vec2<Real>	operator* (const Vec2< Real > & vector,  Real factor ){return Vec2<Real>( factor*vector[0], factor*vector[1] );}
template<class Real>
inline const Vec2<Real>	operator- ( const Vec2<Real> & vector ){ return Vec2< Real>( -vector[0], -vector[1] ); }
template<class Real>
inline const  Vec2<Real> operator- ( const  Vec2<Real> & v1, const  Vec2<Real> & v2 ){ return Vec2< Real >(v1[0]-v2[0], v1[1]-v2[1]); }
template<class Real>
inline const Vec2<Real>	operator/ ( const Vec2<Real> & vector, Real divisor ){ return Vec2<Real>( vector[0]/divisor, vector[1]/divisor ); }
template<class Real> inline bool operator== ( const Vec2< Real >& v1, const Vec2< Real >& v2 ){ return ( ( v1[0]==v2[0] ) && ( v1[1]==v2[1]) ); }
template<class Real> inline bool operator!= ( const Vec2< Real >& v1, const Vec2< Real >& v2 ){ return ( ( v1[0]!=v2[0] ) && ( v1[1]!=v2[1]) ); }

#endif // VEC_H

