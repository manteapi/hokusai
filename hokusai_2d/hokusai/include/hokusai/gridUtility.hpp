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

        void get9Neighbors(std::vector<int>& neighbors, const Vec2d& p, const double radius);
        void get9Neighbors(std::vector<Vec2i>& neighbors,const Vec2d& p, const double radius);

        void get9Neighbors(std::vector<int>& neighbors, const Vec2i& p, const int radius);
        void get9Neighbors(std::vector<Vec2i>& neighbors, const Vec2i& p, const int radius);


        void get5Neighbors(std::vector<int>& neighbors, const Vec2d& p) const;
        void get5Neighbors(std::vector<Vec2i>& neighbors,const Vec2d& p) const;

        void get5Neighbors(std::vector<int>& neighbors, const Vec2i& p) const;
        void get5Neighbors(std::vector<Vec2i>& neighbors, const Vec2i& p) const;

        bool isInside(int id) const;
        bool isInside(int i, int j) const;
        bool isInside(const Vec2i& v) const;
        bool isInside(const Vec2d& v) const;

        Vec2d gridToWorld(const Vec2i& v) const;
        Vec2d gridToWorld(const int i) const;

        Vec2i worldToGrid(const Vec2d& v) const;

        int cellId(const Vec2d& v) const;
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
        double spacing() const;

        void update(const Vec2d& _offset, const Vec2i& _scale, const double& _spacing);
        void update(const Vec2d& _offset, const Vec2d& _scale, const double& _spacing);

        void init(const Vec2d& _offset, const Vec2d& _scale, const double& _spacing);
        void init(const Vec2d& _offset, const Vec2i& _scale, const double& _spacing);
        Vec2i gridCoord(int i) const;
        void info();
        double length(int cellId1, int cellId2) const;

    /*!
     * \brief Return the world coordinates of the four pixel's corners.
     *
     * \return Coordinates of the four pixel's corners.
     */
    std::array< Vec2d, 4 > corners(int pixelId) const;
};

}

#endif
