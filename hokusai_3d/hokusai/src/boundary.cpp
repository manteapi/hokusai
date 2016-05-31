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

#include "../include/hokusai/boundary.hpp"

namespace hokusai
{

Boundary::Boundary()
{
    psi=0.0;
    x=Vec3r(0.0);
    v=Vec3r(0.0);
}

Boundary::Boundary(const Vec3r& _x, const Vec3r _v, const HReal _psi, const BoundaryParams &boundaryParams)
{
    x=_x;
    v=_v;
    psi=_psi;
    m_boundaryParams = boundaryParams;
}

Boundary::Boundary(const Boundary& b)
{
    x=b.x;
    v=b.v;
    psi=b.psi;
    m_boundaryParams = b.boundaryParams();
}

Boundary::~Boundary(){}

BoundaryParams& Boundary::boundaryParams()
{
    return m_boundaryParams;
}

const BoundaryParams& Boundary::boundaryParams() const
{
    return m_boundaryParams;
}

}//namespace hokusai
