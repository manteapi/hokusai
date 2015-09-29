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

#include "../include/hokusai/utils.hpp"

namespace hokusai
{

AkinciKernel::AkinciKernel()
{
    h = 0;
    m_v1 = 0;
    m_v2 = 0;
}

AkinciKernel::AkinciKernel( double _h)
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

MonaghanKernel::MonaghanKernel()
{
    h = 0;
    m_v = 0;
    m_g = 0;
}

MonaghanKernel::MonaghanKernel( double _h )
{
    h = _h;
    m_v = 15.0/(14.0*M_PI*h*h);
    m_g = 15.0/(14.0*M_PI*h*h);
}

MonaghanKernel::MonaghanKernel( const MonaghanKernel& k)
{
    h = k.h;
    m_v = k.m_v;
    m_g = k.m_g;
}

MonaghanKernel::~MonaghanKernel(){}


BoundaryKernel::BoundaryKernel()
{
    h = 0;
    cs = 0;
}

BoundaryKernel::~BoundaryKernel(){}

BoundaryKernel::BoundaryKernel( double _h, double _cs )
{
    h = _h;
    cs = _cs;
}

}
