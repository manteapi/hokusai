#ifndef ALJABR_VEC_HPP
#define ALJABR_VEC_HPP

#include <cmath>
#include <iostream>
#include <array>
#include <algorithm>
#include "Utils.hpp"

namespace aljabr
{
/*
template<class SReal>
class Vec3
{
    private:
        std::array< SReal, 3 > vector;

    public:
        inline Vec3(){ vector.fill(0);}
        inline Vec3(SReal x){vector.fill(x);}
        inline Vec3(const Vec3& v){ vector = v.getConstArray();  }
        inline Vec3( SReal x, SReal y, SReal z){ vector[0] = x; vector[1] = y; vector[2] = z; }
        ~Vec3(){}

        inline void fill(SReal x){vector.fill(x);}
        inline void info() const{ std::cout << "Vec : " << vector[0] << ", " << vector[1] << ", " << vector[2] << std::endl; }

        inline const std::array<SReal, 3>& getConstArray() const { return vector; }
        inline std::array<SReal, 3>& getArrayValue(){ return vector; }
        inline SReal lengthSquared() const { return (vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);}
        inline SReal length() const { return sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);}
        inline void normalize(){ SReal l = sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]); vector[0]/=l; vector[1]/=l; vector[2]/=l;}

        inline Vec3< SReal >& operator= (const Vec3<SReal> & v){ vector = v.getConstArray(); return *this;}

        inline Vec3< SReal > operator-(const Vec3<SReal> & v) const{ return Vec3<SReal>( vector[0]-v[0], vector[1]-v[1], vector[2]-v[2] ); }
        inline Vec3< SReal > operator+(const Vec3<SReal> & v) const{ return Vec3<SReal>( vector[0]+v[0], vector[1]+v[1], vector[2]+v[2] ); }

        inline Vec3< SReal > operator*=(SReal factor) { vector[0]*=factor; vector[1]*=factor; vector[2]*=factor; return *this; }
        inline Vec3< SReal > operator/=(SReal factor) { vector[0]/=factor; vector[1]/=factor; vector[2]/=factor; return *this; }

        inline Vec3< SReal > operator+=(const Vec3<SReal> & v) { vector[0]+=v[0]; vector[1]+=v[1]; vector[2]+=v[2]; return *this; }
        inline Vec3< SReal > operator-=(const Vec3<SReal> & v) { vector[0]-=v[0]; vector[1]-=v[1]; vector[2]-=v[2]; return *this; }

        inline SReal& operator[](int i){ return vector[i]; }
        inline const SReal& operator[](int i) const{ return vector[i]; }

        inline static SReal dotProduct ( const Vec3<SReal>& v1, const Vec3<SReal>& v2){return SReal(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] );}
        inline static Vec3< SReal> crossProduct(const Vec3<SReal>& v1, const Vec3<SReal>& v2) {return Vec3<SReal>(v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]);}

        inline Vec3< SReal > cwiseProd(const Vec3<SReal> &v) const { return Vec3<SReal>( vector[0]*v[0], vector[1]*v[1], vector[2]*v[2]);}

        inline Vec3< SReal > max( SReal value ){ return Vec3<SReal>( std::max( vector[0], value ), std::max( vector[1], value ), std::max( vector[2], value ));}

        inline Vec3< SReal > min( SReal value ){ return Vec3<SReal>( std::min( vector[0], value ), std::min( vector[1], value ), std::min( vector[2], value ));}

        inline static Vec3< SReal > abs( const Vec3<SReal>& v){return Vec3<SReal>(std::abs(v[0]), std::abs(v[1]), std::abs(v[2]));}
        inline static Vec3< SReal > sign( const Vec3<SReal>& v){return Vec3<SReal>(signum<SReal>(v[0]), signum<SReal>(v[1]), signum<SReal>(v[2]));}

        inline static SReal max(Vec3<SReal> v){ return (*std::max_element(v.getConstArray().begin(), v.getConstArray().end())); }
        inline static SReal min(Vec3<SReal> v){ return (*std::min_element(v.getConstArray().begin(), v.getConstArray().end())); }
};

template<class SReal>
inline const Vec3<SReal>	operator* ( SReal factor, const Vec3< SReal > & vector ){return Vec3<SReal>( factor*vector[0], factor*vector[1], factor*vector[2] );}
template<class SReal>
inline const Vec3<SReal>	operator* (const Vec3< SReal > & vector,  SReal factor ){return Vec3<SReal>( factor*vector[0], factor*vector[1], factor*vector[2] );}
template<class SReal>
inline const Vec3<SReal>	operator- ( const Vec3<SReal> & vector ){ return Vec3< SReal>( -vector[0], -vector[1],-vector[2]); }
template<class SReal>
inline const  Vec3<SReal> operator- ( const  Vec3<SReal> & v1, const  Vec3<SReal> & v2 ){ return Vec3< SReal >(v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]); }
template<class SReal>
inline const Vec3<SReal>	operator/ ( const Vec3<SReal> & vector, SReal divisor ){ return Vec3<SReal>( vector[0]/divisor, vector[1]/divisor, vector[2]/divisor ); }
template<class SReal> inline bool operator== ( const Vec3< SReal >& v1, const Vec3< SReal >& v2 ){ return ( ( v1[0]==v2[0] ) && ( v1[1]==v2[1]) && ( v1[2]==v2[2]) ); }
template<class SReal> inline bool operator!= ( const Vec3< SReal >& v1, const Vec3< SReal >& v2 ){ return ( ( v1[0]!=v2[0] ) && ( v1[1]!=v2[1]) && ( v1[2]!=v2[2]) ); }

// Read from an input stream
template<class SReal>
std::istream& operator >> ( std::istream& in, Vec3<SReal>& v )
{
    for( int i=0; i<3; ++i )
    {
        in>>v[i];
    }
    return in;
}

/// Write to an output stream
template<class SReal>
std::ostream& operator << ( std::ostream& out, const Vec3<SReal>& v )
{
for( int i=0; i<3-1; ++i )
{
    out<<v[i]<<" ";
}
out<<v[3-1];
return out;
}


//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//						Vec2< SReal >
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

template<class SReal>
class Vec2
{
    private:
        std::array<SReal, 2> vector;

    public:
        inline Vec2(){vector.fill(0);}
        inline Vec2(SReal x){vector.fill(x);}
        inline Vec2(const Vec2& v){ vector = v.getConstArray();  }
        inline Vec2( SReal x, SReal y){ vector[0] = x; vector[1] = y; }
        ~Vec2(){}

        inline void fill(SReal x){vector.fill(x);}
        inline void info() const{ std::cout << "Vec : " << vector[0] << ", " << vector[1] << std::endl; }

        inline const std::array<SReal, 2>& getConstArray() const { return vector; }
        inline std::array<SReal, 2>& getArrayValue(){ return vector; }
        inline SReal lengthSquared() const { return (vector[0]*vector[0]+vector[1]*vector[1]);}
        inline SReal length() const { return sqrt(vector[0]*vector[0]+vector[1]*vector[1]);}
        inline void normalize(){ SReal l = sqrt(vector[0]*vector[0]+vector[1]*vector[1]); vector[0]/=l; vector[1]/=l;}

        inline Vec2< SReal >& operator= (const Vec2<SReal> & v){ vector = v.getConstArray(); return *this;}

        inline Vec2< SReal > operator-(const Vec2<SReal> & v) const{ return Vec2<SReal>( vector[0]-v[0], vector[1]-v[1] ); }
        inline Vec2< SReal > operator+(const Vec2<SReal> & v) const{ return Vec2<SReal>( vector[0]+v[0], vector[1]+v[1] ); }

        inline Vec2< SReal > operator*=(SReal factor) { vector[0]*=factor; vector[1]*=factor; return *this; }
        inline Vec2< SReal > operator/=(SReal factor) { vector[0]/=factor; vector[1]/=factor; return *this; }

        inline Vec2< SReal > operator+=(const Vec2<SReal> & v) { vector[0]+=v[0]; vector[1]+=v[1]; return *this; }
        inline Vec2< SReal > operator-=(const Vec2<SReal> & v) { vector[0]-=v[0]; vector[1]-=v[1]; return *this; }

        //Copied from Sofa fixed_array.h
        inline bool operator < (const Vec2<SReal>& v ) const
        {
            for( int i=0; i<2; i++ )
            {   
                if( vector[i]<v[i] )
                    return true;  // (*this)<v
                else if( vector[i]>v[i] )
                    return false; // (*this)>v
            }   
            return false; // (*this)==v
        }   


        inline SReal& operator[](int i){ return vector[i]; }
        inline const SReal& operator[](int i) const{ return vector[i]; }

        inline static SReal dotProduct ( const Vec2<SReal>& v1, const Vec2<SReal>& v2){return( v1[0]*v2[0]+v1[1]*v2[1]);}
        inline static SReal crossProduct(const Vec2<SReal>& v1, const Vec2<SReal>& v2){return (v1[0]*v2[1] - v1[1]*v2[0]);}
        inline Vec2< SReal > cwiseProd(const Vec2<SReal> &v) const { return Vec2<SReal>( vector[0]*v[0], vector[1]*v[1]);}
        inline Vec2< SReal > maxVector( SReal value ){ return Vec2<SReal>( std::max( vector[0], value ), std::max( vector[1], value ));}
        inline Vec2< SReal > minVector( SReal value ){ return Vec2<SReal>( std::min( vector[0], value ), std::min( vector[1], value ));}

        inline static SReal max(Vec2<SReal> v){ return (*std::max_element(v.getConstArray().begin(), v.getConstArray().end())); }
        inline static SReal min(Vec2<SReal> v){ return (*std::min_element(v.getConstArray().begin(), v.getConstArray().end())); }
        inline static Vec2< SReal > abs( const Vec2<SReal>& v){return Vec2<SReal>(std::abs(v[0]), std::abs(v[1]));}
        inline static Vec2< SReal > sign( const Vec2<SReal>& v){return Vec2<SReal>(signum<SReal>(v[0]), signum<SReal>(v[1]));}
};

template<class SReal>
inline const Vec2<SReal>	operator* ( SReal factor, const Vec2< SReal > & vector ){return Vec2<SReal>( factor*vector[0], factor*vector[1] );}
template<class SReal>
inline const Vec2<SReal>	operator* (const Vec2< SReal > & vector,  SReal factor ){return Vec2<SReal>( factor*vector[0], factor*vector[1] );}
template<class SReal>
inline const Vec2<SReal>	operator- ( const Vec2<SReal> & vector ){ return Vec2< SReal>( -vector[0], -vector[1] ); }
template<class SReal>
inline const  Vec2<SReal> operator- ( const  Vec2<SReal> & v1, const  Vec2<SReal> & v2 ){ return Vec2< SReal >(v1[0]-v2[0], v1[1]-v2[1]); }
template<class SReal>
inline const Vec2<SReal>	operator/ ( const Vec2<SReal> & vector, SReal divisor ){ return Vec2<SReal>( vector[0]/divisor, vector[1]/divisor ); }
template<class SReal> inline bool operator== ( const Vec2< SReal >& v1, const Vec2< SReal >& v2 ){ return ( ( v1[0]==v2[0] ) && ( v1[1]==v2[1]) ); }
template<class SReal> inline bool operator!= ( const Vec2< SReal >& v1, const Vec2< SReal >& v2 ){ return ( ( v1[0]!=v2[0] ) && ( v1[1]!=v2[1]) ); }

// Read from an input stream
template<class SReal>
std::istream& operator >> ( std::istream& in, Vec2<SReal>& v )
{
    for( int i=0; i<2; ++i )
    {
        in>>v[i];
    }
    return in;
}

/// Write to an output stream
template<class SReal>
std::ostream& operator << ( std::ostream& out, const Vec2<SReal>& v )
{
for( int i=0; i<2-1; ++i )
{
    out<<v[i]<<" ";
}
out<<v[2-1];
return out;
}
*/
}//aljabr

#endif // ALJABR_VEC_HPP
