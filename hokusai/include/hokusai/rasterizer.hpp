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

/*! \file rasterizer.h
 * \brief Useful functions to perform simple rasterization.
 */
#ifndef HOKUSAI_RASTERIZER_H
#define HOKUSAI_RASTERIZER_H

#include <stdio.h>
#include <vector>
#include "Vec.hpp"

namespace hokusai
{
typedef Vec3<int> Vec3i;

//! Rasterize a line between two integers points. 
/*!
 * \brief Points are supposed to have positive coordinates. The code is a C++ translation from the matlab code of Jimmy Shen : <a href="http://www.mathworks.com/matlabcentral/fileexchange/21057-3d-bresenham-s-line-generation">Jimmy Shen</a>

 * \param p1 first point.
 * \param p2 second point.
 * \return A set of integer points which lie between p1 and p2.
 */
std::vector<Vec3i> rasterizeLine(const Vec3i& p1, const Vec3i& p2);


//! Rasterize the span between two lines along a given axis and stop at the shortest line.
/*!
 * \brief minLine's and maxLine's first point should be the same
 * \param minLine shortest line.
 * \param maxLine longest line.
 * \param axis rasterization axis.
 * \param span resulting span as a set of integer points.
 * \param leftover the two points where rasterization ended.
 */
void rasterizeSpan(const std::vector<Vec3i>& minLine, const std::vector<Vec3i>& maxLine, const int& axis,std::vector<Vec3i>& result, std::array< Vec3i, 2>& leftover);


//! Rasterize a triangle.
/*!
 * \param a first point.
 * \param b second point.
 * \param c third point.
 * \return A set of integer points which lie in the triangle.
 */
std::vector<Vec3i> rasterizeTriangle(const Vec3i& a, const Vec3i& b, const Vec3i &c);

//! Rasterize a triangle without using span rasterization.
/*!
 * \param a first point.
 * \param b second point.
 * \param c third point.
 * \return A set of integer points which lie in the triangle.
 */
std::vector<Vec3i> rasterizeTriangleWithoutSpan(const Vec3i& a, const Vec3i& b, const Vec3i &c);

//! Overasterize a triangle.
/*!
 * \brief Sometimes badly shape triangles can leave some holes. This function tries to remove thoses holes by applying several rasterizations of the triangle but with different orientation of the vertices.
 * \param a first point.
 * \param b second point.
 * \param c third point.
 * \return A set of integer points which lie in the triangle.
 */
std::vector<Vec3i> overRasterizeTriangle(const Vec3i& a, const Vec3i& b, const Vec3i &c);

//! Overasterize a triangle without using span rasterization.
/*!
 * \brief Sometimes badly shape triangles can leave some holes. This function tries to remove thoses holes by applying several rasterizations of the triangle but with different orientation of the vertices.
 * \param a first point.
 * \param b second point.
 * \param c third point.
 * \return A set of integer points which lie in the triangle.
 */
std::vector<Vec3i> overRasterizeTriangleWithoutSpan(const Vec3i& a, const Vec3i& b, const Vec3i &c);

//! Rasterize a set of triangles.
/*!
 * \param triangles a set of triangles to rasterize.
 * \return A set which contains for each triangle the points which lie in it.
 */
std::vector<Vec3i> rasterizeTriangles( const std::vector< std::array<Vec3i,3> >& triangles );

//! Rasterize a set of triangles without using span rasterization.
/*!
 * \param triangles a set of triangles to rasterize.
 * \return A set which contains for each triangle the points which lie in it.
 */
std::vector<Vec3i> rasterizeTrianglesWithoutSpan( const std::vector< std::array<Vec3i,3> >& triangles );

//! Rasterize a set of triangles using overrasterization.
/*!
 * \param triangles a set of triangles to rasterize.
 * \return A set which contains for each triangle the points which lie in it.
 */
std::vector<Vec3i> overRasterizeTriangles( const std::vector< std::array<Vec3i,3> >& triangles );

//! Rasterize a set of triangles using overrasterization but without span rasterization.
/*!
 * \param triangles a set of triangles to rasterize.
 * \return A set which contains for each triangle the points which lie in it.
 */
std::vector<Vec3i> overRasterizeTrianglesWithoutSpan( const std::vector< std::array<Vec3i,3> >& triangles );
}
#endif
