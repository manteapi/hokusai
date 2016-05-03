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

#ifndef HOKUSAI_BOUNDARY_HPP
#define HOKUSAI_BOUNDARY_HPP

#include <vector>
#include "common.hpp"

/*! \brief Brief description.
 *         Brief description continued.
 *
 *  Detailed description starts here.
 */

namespace hokusai
{

    class Boundary
    {
        public :
        HReal psi; //density number
        Vec3r x,v; ///position and velocity

        public :
        Boundary()
        {
            psi=0.0;
            x=Vec3r(0.0);
            v=Vec3r(0.0);
        }

        Boundary(const Vec3r& _x, const Vec3r _v = Vec3r(0,0,0), const HReal _psi=0.0)
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
}
#endif //HOKUSAI_BOUNDARY_HPP
