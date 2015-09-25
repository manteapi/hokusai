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

#include <aljabr/Vec.hpp>
#include "triMesh.hpp"
#include <vector>

namespace hokusai
{

typedef aljabr::Vec3<double> Vec3r;
typedef aljabr::Vec2<double> Vec2r;
typedef aljabr::Vec2<int> Vec2i;

bool LineLineIntersect(const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const Vec3r& p4, Vec3r& pa, Vec3r& pb,double & mua, double & mub);
bool LineIntersect(const Vec2r& p1, const Vec2r& p2, const Vec2r& p3, const Vec2r& p4, Vec2r& p);
bool LineIntersect(const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const Vec3r& p4, Vec3r& p);

//Surface sampling
bool AkinciTriangleSampling( const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const double& particleDiameter, std::vector< Vec3r >& samples);
bool AkinciEdgeSampling( const Vec3r& p1, const Vec3r& p2, const double& particleDiameter, std::vector< Vec3r >& samples);
bool AkinciMeshSampling(const TriMesh& mesh, const double &particleDiameter, std::vector<Vec3f>& samples);
bool AkinciFullTriangleSampling( const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const double& particleDiameter, std::vector< Vec3r >& samples);



std::vector<Vec3r > getDiskSampling(const Vec3r& center, double radius, double spacing);
std::vector<Vec3r > getSphereSampling(const Vec3r& center, double radius, double spacingX, double spacingY);
std::vector<Vec3r > getHemiSphereSampling(const Vec3r& center, double radius, double spacingX, double spacingY);

std::vector<Vec3r > getEllipsoidSampling(const Vec3r& center, double axis_1, double axis_2, double axis_3, double spacingX, double spacingY);
std::vector< Vec3r > getTorusSampling(const Vec3r& center, double tubeRadius, double innerRadius, double spacingX, double spacingY);

std::vector< Vec3r > getConeSampling(const Vec3r& center, double height, double stopHeight, double baseRadius, double spacingX, double spacingY);
std::vector<Vec3r > getPyramidSampling(const Vec3r& _center, double base, double height, double spacing);
std::vector< Vec3r > getCylinderSampling(const Vec3r& center, double height, double baseRadius, double spacingX, double spacingY);
std::vector<Vec3r > getCapsuleSampling(const Vec3r& center, double radius, double height, double spacingX, double spacingY);

}
#endif
