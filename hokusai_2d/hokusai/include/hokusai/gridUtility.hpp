#ifndef GRID_UTILITY_HPP
#define GRID_UTILITY_HPP

#include "utility.hpp"
#include <iostream>
#include <array>

namespace hokusai
{

class Grid2dUtility
{
    public :
        Grid2dUtility();
        Grid2dUtility( const Grid2dUtility& gridInfo );
        Grid2dUtility(const Vec2d& _offset, const Vec2d& _scale, const double& _spacing);
        Grid2dUtility(const Vec2d& _offset, const Vec2i& _dimension, const double& _spacing);
        ~Grid2dUtility();

        Vec2d offset;
        Vec2d scale;
        double h;
        Vec2i dimension;

        inline void get9Neighbors(std::vector<int>& neighbors, const Vec2d& p, const double radius)
        {
            get9Neighbors(neighbors, worldToGrid(p), (int)(std::floor(radius/h)));
        }

        inline void get9Neighbors(std::vector<Vec2i>& neighbors,const Vec2d& p, const double radius)
        {
            get9Neighbors(neighbors, worldToGrid(p), (int)std::floor(radius/h) );
        }

        inline void get9Neighbors(std::vector<Vec2i>& neighbors, const Vec2i& p, const int radius)
        {
            neighbors.clear();
            if(!isInside(p))
                return;

                for(int j = std::max(p[1]-radius,0); j<= std::min(p[1]+radius,dimension[1]-1); ++j )
                {
                    for(int i = std::max(p[0]-radius,0); i<= std::min(p[0]+radius,dimension[0]-1); ++i )
                    {
                        neighbors.push_back(Vec2i(i,j));
                    }
                }
        }

        inline void get9Neighbors(std::vector<int>& neighbors, const Vec2i& p, const int radius)
        {
            neighbors.clear();
            if(!isInside(p))
                return;

                for(int j = std::max(p[1]-radius,0); j<= std::min(p[1]+radius,dimension[1]-1); ++j )
                {
                    for(int i = std::max(p[0]-radius,0); i<= std::min(p[0]+radius,dimension[0]-1); ++i )
                    {
                        neighbors.push_back(cellId(i,j));
                    }
                }
        }

        inline bool isInside(int id) const
        {
            return (id>=0 && id<size());
        }

        inline bool isInside(const Vec2d& v) const
        {
            return isInside(worldToGrid(v));
        }

        inline bool isInside(int i, int j) const
        {
            return (i>=0 && i<dimension[0]) && (j>=0 && j<dimension[1]);
        }

        inline bool isInside(const Vec2i& v) const
        {
            return (v[0]>=0 && v[0]<dimension[0]) && (v[1]>=0 && v[1]<dimension[1]);
        }

        inline Vec2d gridToWorld(const int i) const
        {
            return gridToWorld(gridCoord(i));
        }

        inline Vec2d gridToWorld(const Vec2i& v) const
        {
            Vec2d result;
            for(int i=0; i<2; ++i)
                result[i] = offset[i] + h*v[i];
            return result;
        }

        inline Vec2i worldToGrid(const Vec2d& v) const
        {
            Vec2i result;
            for(int i=0; i<2; ++i)
                result[i] = std::floor((v[i]-offset[i])/h);
            return result;
        }

        inline int cellId(const Vec2d& v) const
        {
            return cellId(worldToGrid(v));
        }

        inline int cellId(const Vec2i& v) const
        {
            return v[0] + v[1]*dimension[0];
        }

        inline int cellId(int i, int j) const
        {
            return i + j*dimension[0];
        }

        inline Vec2i gridCoord(int i) const
        {
            Vec2i coord;
            coord[1] = i/dimension[0];
            coord[0] = (i - coord[1]*dimension[0]);
            return coord;
        }

        inline int width() const { return dimension[0]; }
        inline int height() const { return dimension[1]; }
        inline int size() const { return dimension[0]*dimension[1]; }
        inline double spacing() const { return h; }

        void update(const Vec2d& _offset, const Vec2i& _scale, const double& _spacing);
        void update(const Vec2d& _offset, const Vec2d& _scale, const double& _spacing);

        void init(const Vec2d& _offset, const Vec2d& _scale, const double& _spacing);
        void init(const Vec2d& _offset, const Vec2i& _scale, const double& _spacing);
        void info();
};

}

#endif
