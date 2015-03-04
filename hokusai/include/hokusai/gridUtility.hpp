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

#ifndef HOKUSAI_GRID_UTILITY_HPP
#define HOKUSAI_GRID_UTILITY_HPP

#include "Vec.hpp"

namespace hokusai
{
typedef Vec3<double> Vec3f;
typedef Vec3<int> Vec3i;
class GridUtility
{
    public :
        GridUtility();
        GridUtility(const Vec3f& _offset, const Vec3f& _scale, const float& _spacing);
        GridUtility(const Vec3f& _offset, const Vec3i& _dimension, const float& _spacing);
        ~GridUtility();

        Vec3f offset;
        Vec3f scale;
        float h;
        Vec3i dimension;

        void get27Neighbors(std::vector<int>& neighbors, const Vec3f& p, const float radius);
        void get27Neighbors(std::vector<Vec3i>& neighbors,const Vec3f& p, const float radius);

        void get27Neighbors(std::vector<int>& neighbors, const Vec3i& p, const int radius);
        void get27Neighbors(std::vector<Vec3i>& neighbors, const Vec3i& p, const int radius);


        void get7Neighbors(std::vector<int>& neighbors, const Vec3f& p) const;
        void get7Neighbors(std::vector<Vec3i>& neighbors,const Vec3f& p) const;

        void get7Neighbors(std::vector<int>& neighbors, const Vec3i& p) const;
        void get7Neighbors(std::vector<Vec3i>& neighbors, const Vec3i& p) const;

        bool isInside(int id) const;
        bool isInside(int i, int j, int k) const;
        bool isInside(const Vec3i& v) const;
        bool isInside(const Vec3f& v) const;

        Vec3f gridToWorld(const Vec3i& v) const;
        Vec3f gridToWorld(const int i) const;

        Vec3i worldToGrid(const Vec3f& v) const;

        int cellId(const Vec3f& v) const;
        int cellId(const Vec3i& v) const;
        int cellId(int i, int j, int k) const;

        int width() const;
        int height() const;
        int depth() const;
        int size() const;
        float spacing() const;

        void update(const Vec3f& _offset, const Vec3i& _scale, const float& _spacing);
        void update(const Vec3f& _offset, const Vec3f& _scale, const float& _spacing);

        void init(const Vec3f& _offset, const Vec3f& _scale, const float& _spacing);
        void init(const Vec3f& _offset, const Vec3i& _scale, const float& _spacing);
        Vec3i gridCoord(int i) const;
        void info();
};
}
#endif
