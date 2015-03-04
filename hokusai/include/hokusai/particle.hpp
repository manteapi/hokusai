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

#ifndef HOKUSAI_PARTICLE_H
#define HOKUSAI_PARTICLE_H

#include "Vec.hpp"
#include <vector>

/*! \brief Brief description.
 *         Brief description continued.
 *
 *  Detailed description starts here.
 */

namespace hokusai
{

    class Boundary
    {
        typedef Vec3<double> Vec;

        public :
        double psi; //density number
        Vec x,v; ///position and velocity

        public :
        Boundary()
        {
            psi=0.0;
            x=Vec(0.0);
            v=Vec(0.0);
        }

        Boundary(const Vec& _x, const Vec _v = Vec(0,0,0), const double _psi=0.0)
        {
            x=_x;
            v=_v;
            psi=_psi;
        }

        Boundary(const Boundary& b)
        {
            x=b.x;
            v=b.v;
            psi=b.psi;
        }

        ~Boundary(){}
    };

    class Particle
    {
        typedef Vec3<double> Vec;

        public :
        double rho, rho_adv, rho_corr, p, p_l, previousP, aii;
        Vec x, v, v_adv, f_adv, f_p, dii_fluid, dii_boundary, sum_dij, n;
        Vec3<float> c;
        std::vector<int> fluidNeighbor;
        std::vector<int> boundaryNeighbor;

        public :
        Particle()
        {
            rho = rho_adv = rho_corr = p = p_l = previousP = aii = 0.0;
            x = v = v_adv = f_adv = f_p = dii_fluid = dii_boundary = sum_dij = n = Vec(0.0);
            c = Vec3<float>(0.0);
            fluidNeighbor.clear();
            boundaryNeighbor.clear();
        }

        Particle(const Vec& _x, const Vec& _v = Vec(0,0,0), const Vec3<float> _c = Vec3<float>(0,0,1))
        {
            x = _x;
            v = _v;
            c = _c;

            rho = rho_adv = rho_corr = p = p_l = previousP = aii = 0.0;
            v_adv = f_adv = f_p = dii_fluid = dii_boundary = sum_dij = n = Vec(0.0);
            fluidNeighbor.clear();
            boundaryNeighbor.clear();
        }

        Particle(const Particle& _p)
        {
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

            c = _p.c;

            fluidNeighbor=_p.fluidNeighbor;
            boundaryNeighbor=_p.boundaryNeighbor;
        }

        ~Particle(){}
    };
}
#endif
