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

#include <vector>
#include "write_bmp.hpp"
#include "particleContainer.inl"
#include "common.hpp"
#include <fstream>
#include <sstream>

namespace hokusai
{
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

int mortonNumber( std::array<int,3>& index );

//--------------------------------------------------------------------

void transform( HReal p[3], HReal R[3][3] );
void buildRotationMatrix( HReal xrad, HReal yrad, HReal R[3][3] );

//IO functions
void write(const char * filename, std::vector<Vec3r > data);
void write(const char * filename, std::vector<HReal> data);

template<typename ParticleT>
void write_frame(ParticleContainer<ParticleT>& particles, int step, HReal offset=4.0);

//---------------------------------------------------------------------
}
#endif
