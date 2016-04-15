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

#include "../include/hokusai/utils.inl"
#include <bitset>

namespace hokusai
{

AkinciKernel::AkinciKernel()
{
    h = 0;
    m_v1 = 0;
    m_v2 = 0;
}

AkinciKernel::AkinciKernel( HReal _h)
{
    h = _h;
    m_v1 = 32.0/(M_PI*pow(h,9.0));
    m_v2 = pow(h,6.0)/64.0;
    adhesion = 0.007/pow(h,3.25);
}

AkinciKernel::AkinciKernel( const AkinciKernel& k)
{
    h = k.h;
    m_v1 = k.m_v1;
    m_v2 = k.m_v2;
    adhesion = k.adhesion;
}

AkinciKernel::~AkinciKernel(){}

HReal AkinciKernel::cohesionValue( const HReal r )
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

HReal AkinciKernel::adhesionValue( const HReal r)
{
    HReal value=0;
    if( (2.0*r)>h && (r<=h) )
        value=adhesion*pow(-4.0*r*r/h + 6.0*r -2.0*h,1.0/4.0);
    else
        value=0.0;
    return value;
}

MonaghanKernel::MonaghanKernel()
{
    h = 0;
    m_v = 0;
    m_g = 0;
}

MonaghanKernel::MonaghanKernel( HReal _h )
{
    h = _h;
    m_v = 1.0/(4.0*M_PI*h*h*h);
    m_g = 1.0/(4.0*M_PI*h*h*h);
}

MonaghanKernel::MonaghanKernel( const MonaghanKernel& k)
{
    h = k.h;
    m_v = k.m_v;
    m_g = k.m_g;
}

MonaghanKernel::~MonaghanKernel(){}

HReal MonaghanKernel::monaghanValue( const Vec3r & r )
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

void MonaghanKernel::monaghanGradient( const Vec3r& r, Vec3r& gradient )
{
    HReal dist = r.length();
    HReal q = dist/h;
    gradient.fill(0.0);
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

BoundaryKernel::BoundaryKernel()
{
    h = 0;
    cs = 0;
}

BoundaryKernel::~BoundaryKernel(){}

BoundaryKernel::BoundaryKernel( HReal _h, HReal _cs )
{
    h = _h;
    cs = _cs;
}

HReal BoundaryKernel::gamma( HReal distance )
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

int mortonNumber( std::array<int,3>& index )
{
    std::bitset<64> bitsMorton;
    std::array< std::bitset<64>, 3> bitsIndex;
    for(size_t i=0; i<bitsIndex.size(); ++i) bitsIndex[i] = std::bitset<64>(index[i]);
    for(size_t i=0; i<bitsMorton.size()/3; ++i)
        for(size_t j=0; j<bitsIndex.size(); ++j)
            bitsMorton[3*i+j] = bitsIndex[j][i];
    int result = bitsMorton.to_ullong();
    return result;
}

void write(const char * filename, std::vector<HReal> data)
{
    std::ofstream outputFile;
    outputFile.open(filename);
    outputFile.precision(16);
    for(unsigned int i=0; i <data.size(); ++i)
    {
        outputFile << data[i] <<"\n";
    }
    outputFile.close();
}

void write(const char * filename, std::vector<Vec3r > data)
{
    std::ofstream outputFile;
    outputFile.open(filename);
    outputFile.precision(16);
    for(unsigned int i=0; i <data.size(); ++i)
    {
        outputFile << data[i][0] << " " << data[i][1] <<" " <<  data[i][2] << "\n";
    }
    outputFile.close();
}

void buildRotationMatrix( HReal xrad, HReal yrad, HReal R[3][3] ) {
    HReal Rx[3][3] = { {1,0,0},{0,cos(xrad),-sin(xrad)},{0,sin(xrad),cos(xrad)} };
    HReal Ry[3][3] = { {cos(yrad),0,sin(yrad)}, {0,1,0}, {-sin(yrad),0,cos(yrad)} };
    HReal Rtmp[3][3];
    for( int i=0; i<3; i++ ) {
        for( int j=0; j<3; j++ ) {
            R[i][j] = Rtmp[i][j] = Rx[i][j];
        }
    }

    for( int i=0; i<3; i++ ) {
        for( int j=0; j<3; j++ ) {
            R[i][j] = 0.0;
            for( int k=0; k<3; k++ ) R[i][j] += Rtmp[i][k]*Ry[k][j];
        }
    }
}

void transform( HReal p[3], HReal R[3][3] ) {
    HReal p0[3] = { p[0], p[1], p[2] };
    for( int i=0; i<3; i++ ) {
        p[i] = 0.0;
        for( int k=0; k<3; k++ ) {
            p[i] += R[i][k]*p0[k];
        }
    }
}

}//namespace hokusai
