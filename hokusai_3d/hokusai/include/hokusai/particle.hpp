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

#include <vector>
#include "common.hpp"

/*! \brief Brief description.
 *         Brief description continued.
 *
 *  Detailed description starts here.
 */

namespace hokusai
{
    class Particle
    {
        public :
        bool isSurface;
        HReal m, rho, rho_adv, rho_corr, p, p_l, previousP, aii;
        Vec3r x, v, v_adv, f_adv, f_p, dii_fluid, dii_boundary, sum_dij, n;
        std::vector<int> fluidNeighbor;
        std::vector<int> boundaryNeighbor;

        public :
        ~Particle();
        Particle();
        Particle(const Vec3r& _x, const Vec3r& _v);
        Particle(const Particle& _p);

        const HReal& mass() const;
        HReal& mass();
    };
}
#endif
