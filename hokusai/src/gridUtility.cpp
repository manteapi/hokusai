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

#ifndef HOKUSAI_GRID_UTILITY_CPP
#define HOKUSAI_GRID_UTILITY_CPP

#include "../include/hokusai/gridUtility.hpp"

namespace hokusai
{

GridUtility::GridUtility()
{
    offset = Vec3r(0,0,0);
    scale = Vec3r(0,0,0);
    dimension = Vec3i(0,0,0);
    h = 0.0;
}

int GridUtility::width() const { return dimension[0]; }
int GridUtility::height() const { return dimension[1]; }
int GridUtility::depth() const { return dimension[2]; }
int GridUtility::size() const { return dimension[0]*dimension[1]*dimension[2]; }
HReal GridUtility::spacing() const { return h; }

GridUtility::GridUtility(const Vec3r& _offset, const Vec3r& _scale, const HReal& _spacing)
{
    offset = _offset;
    scale = _scale;
    h = _spacing;
    for(int i=0; i<3; ++i)
    {
        dimension[i] = std::floor(scale[i]/h);
    }
}

GridUtility::GridUtility(const Vec3r& _offset, const Vec3i& _dimension, const HReal& _spacing)
{
    offset = _offset;
    dimension = _dimension;
    h = _spacing;
    for(int i=0; i<3; ++i)
    {
        scale[i] = dimension[i]*h;
    }
}

void GridUtility::update(const Vec3r& _offset, const Vec3i& _dimension, const HReal& _spacing)
{
    h = _spacing;
    for(int i=0; i<3; ++i)
    {
        if(_offset[i]<offset[i])
            offset[i] = _offset[i];
        if(_dimension[i]>dimension[i])
        {
            dimension[i] = _dimension[i];
            scale[i] = dimension[i]*h;
        }
    }
}

void GridUtility::update(const Vec3r& _offset, const Vec3r& _scale, const HReal& _spacing)
{
    h = _spacing;
    for(int i=0; i<3; ++i)
    {
        if(_offset[i]<offset[i])
            offset[i] = _offset[i];
        if(_scale[i]>scale[i])
        {
            scale[i] = _scale[i];
            dimension[i] = std::floor(scale[i]/h);
        }
    }
}

void GridUtility::init(const Vec3r& _offset, const Vec3r &_scale, const HReal& _spacing)
{
    offset = _offset;
    scale = _scale;
    h = _spacing;
    for(int i=0; i<3; ++i)
    {
        dimension[i] = std::floor(scale[i]/h);
    }
}

void GridUtility::init(const Vec3r& _offset, const Vec3i& _dimension, const HReal& _spacing)
{
    offset = _offset;
    dimension = _dimension;
    h = _spacing;
    for(int i=0; i<3; ++i)
    {
        scale[i] = dimension[i]*h;
    }
}


GridUtility::~GridUtility(){}

bool GridUtility::isInside(int id) const
{
    if(id>=0 && id<size())
        return true;
    else
        return false;
}


bool GridUtility::isInside(const Vec3r& v) const
{
    Vec3i gridCoord = worldToGrid(v);
    return isInside(gridCoord);
}

bool GridUtility::isInside(int i, int j, int k) const
{
    if( (i>=0 && i<dimension[0]) && (j>=0 && j<dimension[1]) && (k>=0 && k<dimension[2]) )
        return true;
    else
        return false;
}

bool GridUtility::isInside(const Vec3i& v) const
{
    if( (v[0]>=0 && v[0]<dimension[0]) && (v[1]>=0 && v[1]<dimension[1]) && (v[2]>=0 && v[2]<dimension[2]) )
        return true;
    else
        return false;
}

Vec3r GridUtility::gridToWorld(const int i) const
{
    Vec3i gCoord = gridCoord(i);
    Vec3r wCoord = gridToWorld(gCoord);
    return wCoord;
}

Vec3r GridUtility::gridToWorld(const Vec3i& v) const
{
    Vec3r result;
    for(int i=0; i<3; ++i)
        result[i] = offset[i] + h*v[i];
    return result;
}
Vec3i GridUtility::worldToGrid(const Vec3r& v) const
{
    Vec3i result;
    for(int i=0; i<3; ++i)
        result[i] = std::floor((v[i]-offset[i])/h);
    return result;
}
int GridUtility::cellId(const Vec3r& v) const
{
    Vec3i gridCoord = worldToGrid(v);
    return cellId(gridCoord);
}

int GridUtility::cellId(const Vec3i& v) const
{
    return v[0] + v[1]*dimension[0] + v[2]*dimension[0]*dimension[1];
}

int GridUtility::cellId(int i, int j, int k) const
{
    return i + j*dimension[0] + k*dimension[0]*dimension[1];
}

void GridUtility::get27Neighbors(std::vector<int>& neighbors, const int i, const HReal radius) const
{
    Vec3i gCoord= gridCoord(i);
    get27Neighbors(neighbors, gCoord, std::floor(radius/h) );
}

void GridUtility::get27Neighbors(std::vector<Vec3i>& neighbors,const int i, const HReal radius) const
{
    Vec3i gCoord = gridCoord(i);
    get27Neighbors(neighbors, gCoord, std::floor(radius/h) );
}

void GridUtility::get27Neighbors(std::vector<int>& neighbors, const Vec3r& p, const HReal radius) const
{
    Vec3i gridCoord = worldToGrid(p);
    get27Neighbors(neighbors, gridCoord, std::floor(radius/h) );
}

void GridUtility::get27Neighbors(std::vector<Vec3i>& neighbors,const Vec3r& p, const HReal radius) const 
{
    Vec3i gridCoord = worldToGrid(p);
    get27Neighbors(neighbors, gridCoord, std::floor(radius/h) );
}

void GridUtility::get27Neighbors(std::vector<Vec3i>& neighbors, const Vec3i& p, const int radius) const
{
    if(!isInside(p))
    {
        neighbors.clear();
        return;
    }

    int kmin = std::max(p[2]-radius,0);
    int kmax = std::min(p[2]+radius,dimension[2]-1);
    int jmin = std::max(p[1]-radius,0);
    int jmax = std::min(p[1]+radius,dimension[1]-1);
    int imin = std::max(p[0]-radius,0);
    int imax = std::min(p[0]+radius,dimension[0]-1);
    int size = (kmax-kmin+1)*(jmax-jmin+1)*(imax-imin+1);
    neighbors.resize(size);
    int counter = 0;

    for(int k = kmin; k<= kmax; ++k )
    {
        for(int j = jmin; j<= jmax; ++j )
        {
            for(int i = imin; i<= imax; ++i )
            {
                neighbors[counter] = Vec3i(i,j,k);
                counter++;
            }
        }
    }
}

void GridUtility::get27Neighbors(std::vector<int>& neighbors, const Vec3i& p, const int radius) const
{
    if(!isInside(p))
    {
        neighbors.clear();
        return;
    }

    int kmin = std::max(p[2]-radius,0);
    int kmax = std::min(p[2]+radius,dimension[2]-1);
    int jmin = std::max(p[1]-radius,0);
    int jmax = std::min(p[1]+radius,dimension[1]-1);
    int imin = std::max(p[0]-radius,0);
    int imax = std::min(p[0]+radius,dimension[0]-1);
    int size = (kmax-kmin+1)*(jmax-jmin+1)*(imax-imin+1);
    neighbors.resize(size);
    int counter = 0;

    for(int k = kmin; k<= kmax; ++k )
    {
        for(int j = jmin; j<= jmax; ++j )
        {
            for(int i = imin; i<= imax; ++i )
            {
                neighbors[counter] = cellId(i,j,k);
                counter++;
            }
        }
    }
}


void GridUtility::get7Neighbors(std::vector<int>& neighbors, const Vec3r& p) const
{
    Vec3i gridCoord = worldToGrid(p);
    get7Neighbors(neighbors, gridCoord);
}

void GridUtility::get7Neighbors(std::vector<Vec3i>& neighbors,const Vec3r& p) const
{
    Vec3i gridCoord = worldToGrid(p);
    get7Neighbors(neighbors, gridCoord);
}

void GridUtility::get7Neighbors(std::vector<Vec3i>& neighbors, const Vec3i& p) const
{
    neighbors.clear();
    neighbors.reserve(7);
    if(!isInside(p))
        return;

    Vec3i neighbor;

    neighbor = p;
    if(isInside(neighbor))
        neighbors.push_back(neighbor);

    neighbor = p + Vec3i(1,0,0);
    if(isInside(neighbor))
        neighbors.push_back(neighbor);

    neighbor = p + Vec3i(-1,0,0);
    if(isInside(neighbor))
        neighbors.push_back(neighbor);

    neighbor = p + Vec3i(0,1,0);
    if(isInside(neighbor))
        neighbors.push_back(neighbor);

    neighbor = p + Vec3i(0,-1,0);
    if(isInside(neighbor))
        neighbors.push_back(neighbor);

    neighbor = p + Vec3i(0,0,1);
    if(isInside(neighbor))
        neighbors.push_back(neighbor);

    neighbor = p + Vec3i(0,0,-1);
    if(isInside(neighbor))
        neighbors.push_back(neighbor);
}

void GridUtility::get7Neighbors(std::vector<int>& neighbors, const Vec3i& p) const
{
    neighbors.clear();
    neighbors.reserve(7);
    if(!isInside(p))
        return;

    Vec3i neighbor;

    neighbor = p;
    if(isInside(neighbor))
        neighbors.push_back(cellId(neighbor));

    neighbor = p + Vec3i(1,0,0);
    if(isInside(neighbor))
        neighbors.push_back(cellId(neighbor));

    neighbor = p + Vec3i(-1,0,0);
    if(isInside(neighbor))
        neighbors.push_back(cellId(neighbor));

    neighbor = p + Vec3i(0,1,0);
    if(isInside(neighbor))
        neighbors.push_back(cellId(neighbor));

    neighbor = p + Vec3i(0,-1,0);
    if(isInside(neighbor))
        neighbors.push_back(cellId(neighbor));

    neighbor = p + Vec3i(0,0,1);
    if(isInside(neighbor))
        neighbors.push_back(cellId(neighbor));

    neighbor = p + Vec3i(0,0,-1);
    if(isInside(neighbor))
        neighbors.push_back(cellId(neighbor));
}

Vec3i GridUtility::gridCoord(int i) const
{
    Vec3i coord;
    coord[2] = i/(dimension[0]*dimension[1]);
    coord[1] = (i - coord[2]*dimension[0]*dimension[1])/dimension[0];
    coord[0] = (i - coord[1]*dimension[0] - coord[2]*dimension[0]*dimension[1]);
    return coord;
}

void GridUtility::info()
{
    std::cout << "Dimension: " << width() << "x" << height() << "x" << depth() << " = " << size() << std::endl;
    std::cout << "Spacing : " << spacing() << std::endl;
    std::cout << "Offset : " << offset[0] << ", " << offset[1] << ", " << offset[2] << std::endl;
    std::cout << "Scale : " << scale[0] << ", " << scale[1] << ", " << scale[2] << std::endl;
    std::cout << "Max : " << offset[0]+scale[0] << ", " << offset[1]+scale[1] << ", " << offset[2] + scale[2] << std::endl;
}

}
#endif
