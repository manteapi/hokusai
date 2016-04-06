/*
* Copyright 2015 Pierre-Luc Manteaux
*
*This file is part of Hokusai.
*
*Hokusai is free software: you can redistribute it and/or modify
*it under the terms of the GNU General Public License as published by
*the Free Software Foundation, either version 3 of the License, or
*(at your option) any later version.
*
*Hokusai is distributed in the hope that it will be useful,
*but WITHOUT ANY WARRANTY; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*GNU General Public License for more details.
*
*You should have received a copy of the GNU General Public License
*along with Hokusai.  If not, see <http://www.gnu.org/licenses/>.
*
* Author : Pierre-Luc Manteaux
* Contact : pierre-luc.manteaux@inria.fr
*/

#ifndef HOKUSAI_UTILS_HPP
#define HOKUSAI_UTILS_HPP

//#include <GL/gl.h>
#include <vector>
#include <magnet/math/morton_number.hpp>
#include "write_bmp.hpp"
#include "particle.hpp"
#include "common.hpp"
#include <fstream>
#include <sstream>

namespace hokusai
{
using namespace std;

class AkinciKernel
{
public :
    AkinciKernel();
    AkinciKernel(HReal _h);
    AkinciKernel(const AkinciKernel& k);
    ~AkinciKernel();
    HReal cohesionValue( const HReal r);
    HReal adhesionValue(const HReal r);
    HReal h,m_v1,m_v2,adhesion;
};

//Monaghan 3D kernel
class MonaghanKernel
{
public :

    MonaghanKernel();
    MonaghanKernel( HReal _h );
    MonaghanKernel( const MonaghanKernel& k );
    ~MonaghanKernel();

public :

    HReal h;
    HReal m_v;
    HReal m_g;

public :

    HReal monaghanValue( const Vec3r& r );
    void monaghanGradient( const Vec3r& r, Vec3r& gradient );
};

class BoundaryKernel
{
public :

    BoundaryKernel();
    BoundaryKernel( HReal h, HReal cs );
    BoundaryKernel( const BoundaryKernel& k);
    ~BoundaryKernel();

public :

    HReal h;
    HReal cs;

public :

    HReal gamma( HReal distance );
};

inline int mortonNumber( array<int,3>& index )
{
    size_t x = index[0];
    size_t y = index[1];
    size_t z = index[0];
    magnet::math::MortonNumber<3> m(x,y,z);
    size_t morton_num = m.getMortonNum();
    return morton_num;
}

void insertionSort( vector< pair<int,int> >& data );

class Box
{
public :
    Vec3r min;
    Vec3r max;
public :
    Box();
    Box( Vec3r& _min, Vec3r& _max);
    ~Box();
public :
    //void draw();
};

/*
class SolidSphere
{
    protected :
        std::vector<GLHReal> vertices;
        std::vector<GLHReal> normals;
        std::vector<GLHReal> texcoords;
        std::vector<GLuint> indices;

    public :
        ~SolidSphere();
        SolidSphere(HReal radius, unsigned int rings, unsigned int sectors);
    public :
        //void draw(GLHReal x, GLHReal y, GLHReal z);
};
*/

//------------------------------INLINE DEFINITION-----------------------------------
inline HReal MonaghanKernel::monaghanValue( const Vec3r & r )
{
    HReal value = 0.0;
    HReal q = r.length()/h;
    if( q >= 0 && q < 1 )
    {
        value = m_v*( (2-q)*(2-q)*(2-q) - 4.0f*(1-q)*(1-q)*(1-q));
    }
    else if ( q >=1 && q < 2 )
    {
        value = m_v*( (2-q)*(2-q)*(2-q) );
    }
    else
    {
        value = 0.0f;
    }
    return value;
}

inline void MonaghanKernel::monaghanGradient( const Vec3r& r, Vec3r& gradient )
{
    HReal dist = r.length();
    HReal q = dist/h;
    gradient.setAllValue(0.0);
    if( q >= 0 && q < 1 )
    {
        HReal scalar = -3.0f*(2-q)*(2-q);
        scalar += 12.0f*(1-q)*(1-q);
        gradient = (m_g*scalar/(dist*h))*r;
    }
    else if ( q >=1 && q < 2 )
    {
        HReal scalar = -3.0f*(2-q)*(2-q);
        gradient = (m_g*scalar/(dist*h))*r;
    }
}

inline HReal BoundaryKernel::gamma( HReal distance )
{
    HReal q = distance / h;
    HReal coeff = 0.0;
    if( q > 0 && q <= 0.666666667f ) //2.0/3.0
    {
        coeff = 0.666666667f;
    }
    else if( q > 0.666666667f && q < 1 )
    {
        coeff = q*(2 - 1.5*q);
    }
    else if( q >= 1 && q < 2 )
    {
        coeff = 0.5*(2-q)*(2-q);
    }
    coeff *= (0.02*cs*cs/distance);
    return coeff;
}


inline HReal AkinciKernel::cohesionValue( const HReal r )
{
    HReal value=0;
    if( (2.0*r>h) && (r<=h) )
    {
        //value=m_v1*pow(h-r,3.0)*pow(r,3.0);
        value=m_v1*(h-r)*(h-r)*(h-r)*r*r*r;
    }
    else if( (r>0.0) && (2.0*r<=h) )
    {
        //value=m_v1*(2.0*pow(h-r,3.0)*pow(r,3.0)-m_v2);
        value=m_v1*(2.0*(h-r)*(h-r)*(h-r)*r*r*r-m_v2);
    }
    else
        value=0.0;
    return value;
}

inline HReal AkinciKernel::adhesionValue( const HReal r)
{
    HReal value=0;
    if( (2.0*r)>h && (r<=h) )
        value=adhesion*pow(-4.0*r*r/h + 6.0*r -2.0*h,1.0/4.0);
    else
        value=0.0;
    return value;
}

//--------------------------------------------------------------------

//IO functions
void write(const char * filename, vector<Vec3r > data);
void write(const char * filename, vector<HReal> data);
void write_frame(vector<Particle>& particles, int step, HReal offset=4.0);

//---------------------------------------------------------------------
}
#endif
