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
#include <aljabr/Vec.hpp>
#include <magnet/math/morton_number.hpp>
#include "write_bmp.hpp"
#include "particle.hpp"
#include <fstream>
#include <sstream>

namespace hokusai
{
using namespace std;

typedef aljabr::Vec3<double> Vec;

class AkinciKernel
{
public :
    AkinciKernel();
    AkinciKernel(double _h);
    AkinciKernel(const AkinciKernel& k);
    ~AkinciKernel();
    double cohesionValue( const double r);
    double adhesionValue(const double r);
    double h,m_v1,m_v2,adhesion;
};

//Monaghan 3D kernel
class MonaghanKernel
{
    typedef Vec3r Vec;

public :

    MonaghanKernel();
    MonaghanKernel( double _h );
    MonaghanKernel( const MonaghanKernel& k );
    ~MonaghanKernel();

public :

    double h;
    double m_v;
    double m_g;

public :

    double monaghanValue( const Vec& r );
    void monaghanGradient( const Vec& r, Vec& gradient );
};

class BoundaryKernel
{
public :

    BoundaryKernel();
    BoundaryKernel( double h, double cs );
    BoundaryKernel( const BoundaryKernel& k);
    ~BoundaryKernel();

public :

    double h;
    double cs;

public :

    double gamma( double distance );
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
    typedef Vec3r Vec;
public :
    Vec min;
    Vec max;
public :
    Box();
    Box( Vec& _min, Vec& _max);
    ~Box();
public :
    //void draw();
};

/*
class SolidSphere
{
    protected :
        std::vector<GLfloat> vertices;
        std::vector<GLfloat> normals;
        std::vector<GLfloat> texcoords;
        std::vector<GLuint> indices;

    public :
        ~SolidSphere();
        SolidSphere(double radius, unsigned int rings, unsigned int sectors);
    public :
        //void draw(GLfloat x, GLfloat y, GLfloat z);
};
*/

//------------------------------INLINE DEFINITION-----------------------------------
inline double MonaghanKernel::monaghanValue( const Vec & r )
{
    double value = 0.0;
    double q = r.length()/h;
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

inline void MonaghanKernel::monaghanGradient( const Vec& r, Vec& gradient )
{
    double dist = r.length();
    double q = dist/h;
    gradient.setAllValue(0.0);
    if( q >= 0 && q < 1 )
    {
        double scalar = -3.0f*(2-q)*(2-q);
        scalar += 12.0f*(1-q)*(1-q);
        gradient = (m_g*scalar/(dist*h))*r;
    }
    else if ( q >=1 && q < 2 )
    {
        double scalar = -3.0f*(2-q)*(2-q);
        gradient = (m_g*scalar/(dist*h))*r;
    }
}

inline double BoundaryKernel::gamma( double distance )
{
    double q = distance / h;
    double coeff = 0.0;
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


inline double AkinciKernel::cohesionValue( const double r )
{
    double value=0;
    if( (2.0*r>h) && (r<=h) )
        value=m_v1*pow(h-r,3.0)*pow(r,3.0);
    else if( (r>0.0) && (2.0*r<=h) )
        value=m_v1*(2.0*pow(h-r,3.0)*pow(r,3.0)-m_v2);
    else
        value=0.0;
    return value;
}

inline double AkinciKernel::adhesionValue( const double r)
{
    double value=0;
    if( (2.0*r)>h && (r<=h) )
        value=adhesion*pow(-4.0*r*r/h + 6.0*r -2.0*h,1.0/4.0);
    else
        value=0.0;
    return value;
}

//--------------------------------------------------------------------

//IO functions
void write(const char * filename, vector<Vec3r > data);
void write(const char * filename, vector<double> data);
void write_frame(vector<Particle>& particles, int step, float offset=4.0);

//---------------------------------------------------------------------
}
#endif
