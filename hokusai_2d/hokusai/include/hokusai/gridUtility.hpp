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

#include "utility.hpp"

namespace hokusai
{

class Grid2dUtility
{
    public :
        Grid2dUtility();
        Grid2dUtility( const Grid2dUtility& gridInfo );
        Grid2dUtility(const Vec2r& _offset, const Vec2r& _scale, const SReal& _spacing);
        Grid2dUtility(const Vec2r& _offset, const Vec2i& _dimension, const SReal& _spacing);
        ~Grid2dUtility();

        Vec2r offset;
        Vec2r scale;
        SReal h;
        Vec2i dimension;

        void get9Neighbors(std::vector<int>& neighbors, const Vec2r& p, const SReal radius);
        void get9Neighbors(std::vector<Vec2i>& neighbors,const Vec2r& p, const SReal radius);

        void get9Neighbors(std::vector<int>& neighbors, const Vec2i& p, const int radius);
        void get9Neighbors(std::vector<Vec2i>& neighbors, const Vec2i& p, const int radius);


        void get5Neighbors(std::vector<int>& neighbors, const Vec2r& p) const;
        void get5Neighbors(std::vector<Vec2i>& neighbors,const Vec2r& p) const;

        void get5Neighbors(std::vector<int>& neighbors, const Vec2i& p) const;
        void get5Neighbors(std::vector<Vec2i>& neighbors, const Vec2i& p) const;

        bool isInside(int id) const;
        bool isInside(int i, int j) const;
        bool isInside(const Vec2i& v) const;
        bool isInside(const Vec2r& v) const;

        Vec2r gridToWorld(const Vec2i& v) const;
        Vec2r gridToWorld(const int i) const;

        Vec2i worldToGrid(const Vec2r& v) const;

        int cellId(const Vec2r& v) const;
        int cellId(const Vec2i& v) const;
        int cellId(int i, int j) const;

        /*! \brief Return the reference to the neighbor of a pixel givent the side of the neighbor.
        *
        * For the side of the neighbor, 0 corresponds to the bottom side, 1 to the right, 2 to the top and 3 to the left.
        *
        * \param pixelId a reference to the current pixel.
        * \param neighborId the side of the neighbor
        * \return a reference to the neighbor pixel
        */
        int neighborPixelId(int pixelId, int neighborId) const;

        int width() const;
        int height() const;
        int size() const;
        SReal spacing() const;

        void update(const Vec2r& _offset, const Vec2i& _scale, const SReal& _spacing);
        void update(const Vec2r& _offset, const Vec2r& _scale, const SReal& _spacing);

        void init(const Vec2r& _offset, const Vec2r& _scale, const SReal& _spacing);
        void init(const Vec2r& _offset, const Vec2i& _scale, const SReal& _spacing);
        Vec2i gridCoord(int i) const;
        void info();
};

}
#endif
