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

#include "../include/hokusai/particle.hpp"

namespace hokusai
{

Particle::Particle()
{
    isSurface=true;
    rho = rho_adv = rho_corr = p = p_l = previousP = aii = 0.0;
    m = 1.0;
    x = v = v_adv = f_adv = f_p = dii_fluid = dii_boundary = sum_dij = n = Vec3r(0.0);
    fluidNeighbor.clear();
    boundaryNeighbor.clear();
}

Particle::Particle(const Vec3r& _x, const Vec3r& _v, const HReal &_m)
{
    x = _x;
    v = _v;
    m = _m;
    isSurface=true;
    rho = rho_adv = rho_corr = p = p_l = previousP = aii = 0.0;
    v_adv = f_adv = f_p = dii_fluid = dii_boundary = sum_dij = n = Vec3r(0.0);
    fluidNeighbor.clear();
    boundaryNeighbor.clear();
}

Particle::Particle(const Particle& _p)
{
    m = _p.mass();
    rho=_p.rho;
    rho_adv=_p.rho_adv;
    rho_corr=_p.rho_corr;
    p=_p.p;
    p_l=_p.p_l;
    previousP=_p.previousP;
    aii=_p.aii;

    x = _p.x;
    v = _p.v;
    v_adv = _p.v_adv;
    f_adv = _p.f_adv;
    f_p = _p.f_p;
    dii_fluid=_p.dii_fluid;
    dii_boundary=_p.dii_boundary;
    sum_dij=_p.sum_dij;
    n = _p.n;

    isSurface=_p.isSurface;

    fluidNeighbor=_p.fluidNeighbor;
    boundaryNeighbor=_p.boundaryNeighbor;
}

Particle::~Particle(){}


const HReal& Particle::mass() const
{
    return m;
}

HReal& Particle::mass()
{
    m;
}

}//namespace hokusai
