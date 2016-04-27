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
#include "particle.hpp"
#include "common.hpp"
#include <fstream>
#include <sstream>

namespace hokusai
{

int mortonNumber( std::array<int,3>& index );

//--------------------------------------------------------------------

//IO functions
void write(const char * filename, std::vector<Vec3r > data);
void write(const char * filename, std::vector<HReal> data);
void write_frame(std::vector<Particle>& particles, int step, HReal offset=4.0);

//---------------------------------------------------------------------
}
#endif
