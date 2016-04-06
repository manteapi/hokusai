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

#ifndef HOKUSAI_SAMPLER_HPP
#define HOKUSAI_SAMPLER_HPP

#include "triMesh.hpp"
#include "common.hpp"
#include <vector>

namespace hokusai
{

bool LineLineIntersect(const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const Vec3r& p4, Vec3r& pa, Vec3r& pb,HReal & mua, HReal & mub);
bool LineIntersect(const Vec2r& p1, const Vec2r& p2, const Vec2r& p3, const Vec2r& p4, Vec2r& p);
bool LineIntersect(const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const Vec3r& p4, Vec3r& p);

//Surface sampling
bool AkinciTriangleSampling( const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const HReal& particleDiameter, std::vector< Vec3r >& samples);
bool AkinciEdgeSampling( const Vec3r& p1, const Vec3r& p2, const HReal& particleDiameter, std::vector< Vec3r >& samples);
bool AkinciMeshSampling(const TriMesh& mesh, const HReal &particleDiameter, std::vector<Vec3r>& samples);
bool AkinciFullTriangleSampling( const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const HReal& particleDiameter, std::vector< Vec3r >& samples);



std::vector<Vec3r > getDiskSampling(const Vec3r& center, HReal radius, HReal spacing);
std::vector<Vec3r > getSphereSampling(const Vec3r& center, HReal radius, HReal spacingX, HReal spacingY);
std::vector<Vec3r > getHemiSphereSampling(const Vec3r& center, HReal radius, HReal spacingX, HReal spacingY);

std::vector<Vec3r > getEllipsoidSampling(const Vec3r& center, HReal axis_1, HReal axis_2, HReal axis_3, HReal spacingX, HReal spacingY);
std::vector< Vec3r > getTorusSampling(const Vec3r& center, HReal tubeRadius, HReal innerRadius, HReal spacingX, HReal spacingY);

std::vector< Vec3r > getConeSampling(const Vec3r& center, HReal height, HReal stopHeight, HReal baseRadius, HReal spacingX, HReal spacingY);
std::vector<Vec3r > getPyramidSampling(const Vec3r& _center, HReal base, HReal height, HReal spacing);
std::vector< Vec3r > getCylinderSampling(const Vec3r& center, HReal height, HReal baseRadius, HReal spacingX, HReal spacingY);
std::vector<Vec3r > getCapsuleSampling(const Vec3r& center, HReal radius, HReal height, HReal spacingX, HReal spacingY);

}
#endif
