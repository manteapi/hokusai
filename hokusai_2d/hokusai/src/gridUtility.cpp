#ifndef GRID_UTILITY_CPP
#define GRID_UTILITY_CPP

#include "../include/mosaic/gridUtility.hpp"

namespace mosaic
{

Grid2dUtility::Grid2dUtility()
{
    offset = Vec2r(0,0);
    scale = Vec2r(0,0);
    dimension = Vec2i(0,0);
    h = 0.0;
}

Grid2dUtility::Grid2dUtility( const Grid2dUtility& gridInfo )
{
    offset = gridInfo.offset;
    scale = gridInfo.scale;
    dimension = gridInfo.dimension;
    h = gridInfo.h;
}

int Grid2dUtility::width() const { return dimension[0]; }
int Grid2dUtility::height() const { return dimension[1]; }
int Grid2dUtility::size() const { return dimension[0]*dimension[1]; }
SReal Grid2dUtility::spacing() const { return h; }

Grid2dUtility::Grid2dUtility(const Vec2r& _offset, const Vec2r& _scale, const SReal& _spacing)
{
    offset = _offset;
    scale = _scale;
    h = _spacing;
    for(int i=0; i<2; ++i)
    {
        dimension[i] = std::floor(scale[i]/h);
    }
}

Grid2dUtility::Grid2dUtility(const Vec2r& _offset, const Vec2i& _dimension, const SReal& _spacing)
{
    offset = _offset;
    dimension = _dimension;
    h = _spacing;
    for(int i=0; i<2; ++i)
    {
        scale[i] = dimension[i]*h;
    }
}

void Grid2dUtility::update(const Vec2r& _offset, const Vec2i& _dimension, const SReal& _spacing)
{
    h = _spacing;
    for(int i=0; i<2; ++i)
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

void Grid2dUtility::update(const Vec2r& _offset, const Vec2r& _scale, const SReal& _spacing)
{
    h = _spacing;
    for(int i=0; i<2; ++i)
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

void Grid2dUtility::init(const Vec2r& _offset, const Vec2r &_scale, const SReal& _spacing)
{
    offset = _offset;
    scale = _scale;
    h = _spacing;
    for(int i=0; i<2; ++i)
    {
        dimension[i] = std::floor(scale[i]/h);
    }
}

void Grid2dUtility::init(const Vec2r& _offset, const Vec2i& _dimension, const SReal& _spacing)
{
    offset = _offset;
    dimension = _dimension;
    h = _spacing;
    for(int i=0; i<2; ++i)
    {
        scale[i] = dimension[i]*h;
    }
}


Grid2dUtility::~Grid2dUtility(){}

bool Grid2dUtility::isInside(int id) const
{
    if(id>=0 && id<size())
        return true;
    else
        return false;
}


bool Grid2dUtility::isInside(const Vec2r& v) const
{
    Vec2i gridCoord = worldToGrid(v);
    return isInside(gridCoord);
}

bool Grid2dUtility::isInside(int i, int j) const
{
    if( (i>=0 && i<dimension[0]) && (j>=0 && j<dimension[1]) )
        return true;
    else
        return false;
}

bool Grid2dUtility::isInside(const Vec2i& v) const
{
    if( (v[0]>=0 && v[0]<dimension[0]) && (v[1]>=0 && v[1]<dimension[1]) )
        return true;
    else
        return false;
}

Vec2r Grid2dUtility::gridToWorld(const int i) const
{
    Vec2i gCoord = gridCoord(i);
    Vec2r wCoord = gridToWorld(gCoord);
    return wCoord;
}

Vec2r Grid2dUtility::gridToWorld(const Vec2i& v) const
{
    Vec2r result;
    for(int i=0; i<2; ++i)
        result[i] = offset[i] + h*v[i];
    return result;
}
Vec2i Grid2dUtility::worldToGrid(const Vec2r& v) const
{
    Vec2i result;
    for(int i=0; i<2; ++i)
        result[i] = std::floor((v[i]-offset[i])/h);
    return result;
}
int Grid2dUtility::cellId(const Vec2r& v) const
{
    Vec2i gridCoord = worldToGrid(v);
    return cellId(gridCoord);
}

int Grid2dUtility::cellId(const Vec2i& v) const
{
    return v[0] + v[1]*dimension[0];
}

int Grid2dUtility::cellId(int i, int j) const
{
    return i + j*dimension[0];
}

int Grid2dUtility::neighborPixelId(int pixelId, int neighborId) const
{
    Vec2i neighborCoord = gridCoord(pixelId);

    switch( neighborId )
    {
        case 0 :
            neighborCoord += Vec2i(0,-1);
            break;
        case 1 :
            neighborCoord += Vec2i(1,0);
            break;
        case 2 :
            neighborCoord += Vec2i(0,1);
            break;
        case 3 :
            neighborCoord += Vec2i(-1,0);
            break;
        default :
            break;
    }

    if ( isInside(neighborCoord) )
        return cellId( neighborCoord );
    else
        return -1;
}

void Grid2dUtility::get9Neighbors(std::vector<int>& neighbors, const Vec2r& p, const SReal radius)
{
    Vec2i gridCoord = worldToGrid(p);
    get9Neighbors(neighbors, gridCoord, std::floor(radius/h) );
}

void Grid2dUtility::get9Neighbors(std::vector<Vec2i>& neighbors,const Vec2r& p, const SReal radius)
{
    Vec2i gridCoord = worldToGrid(p);
    get9Neighbors(neighbors, gridCoord, std::floor(radius/h) );
}

void Grid2dUtility::get9Neighbors(std::vector<Vec2i>& neighbors, const Vec2i& p, const int radius)
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

void Grid2dUtility::get9Neighbors(std::vector<int>& neighbors, const Vec2i& p, const int radius)
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


void Grid2dUtility::get5Neighbors(std::vector<int>& neighbors, const Vec2r& p) const
{
    Vec2i gridCoord = worldToGrid(p);
    get5Neighbors(neighbors, gridCoord);
}

void Grid2dUtility::get5Neighbors(std::vector<Vec2i>& neighbors,const Vec2r& p) const
{
    Vec2i gridCoord = worldToGrid(p);
    get5Neighbors(neighbors, gridCoord);
}

void Grid2dUtility::get5Neighbors(std::vector<Vec2i>& neighbors, const Vec2i& p) const
{
    neighbors.clear();
    if(!isInside(p))
        return;

    Vec2i neighbor;

    neighbor = p;
    if(isInside(neighbor))
        neighbors.push_back(neighbor);

    neighbor = p + Vec2i(1,0);
    if(isInside(neighbor))
        neighbors.push_back(neighbor);

    neighbor = p + Vec2i(-1,0);
    if(isInside(neighbor))
        neighbors.push_back(neighbor);

    neighbor = p + Vec2i(0,1);
    if(isInside(neighbor))
        neighbors.push_back(neighbor);

    neighbor = p + Vec2i(0,-1);
    if(isInside(neighbor))
        neighbors.push_back(neighbor);
}

void Grid2dUtility::get5Neighbors(std::vector<int>& neighbors, const Vec2i& p) const
{
    neighbors.clear();
    if(!isInside(p))
        return;

    Vec2i neighbor;

    neighbor = p;
    if(isInside(neighbor))
        neighbors.push_back(cellId(neighbor));

    neighbor = p + Vec2i(1,0);
    if(isInside(neighbor))
        neighbors.push_back(cellId(neighbor));

    neighbor = p + Vec2i(-1,0);
    if(isInside(neighbor))
        neighbors.push_back(cellId(neighbor));

    neighbor = p + Vec2i(0,1);
    if(isInside(neighbor))
        neighbors.push_back(cellId(neighbor));

    neighbor = p + Vec2i(0,-1);
    if(isInside(neighbor))
        neighbors.push_back(cellId(neighbor));
}

Vec2i Grid2dUtility::gridCoord(int i) const
{
    Vec2i coord;
    coord[1] = i/dimension[0];
    coord[0] = (i - coord[1]*dimension[0]);
    return coord;
}

void Grid2dUtility::info()
{
    std::cout << "Dimension: " << width() << "x" << height() << "x" << " = " << size() << std::endl;
    std::cout << "Spacing : " << spacing() << std::endl;
    std::cout << "Offset : " << offset[0] << ", " << offset[1] << std::endl;
    std::cout << "Scale : " << scale[0] << ", " << scale[1] << std::endl;
    std::cout << "Max : " << offset[0]+scale[0] << ", " << offset[1]+scale[1] << std::endl;
}

SReal Grid2dUtility::length(int cellId1, int cellId2) const
{
    Vec2r wCoord1 = gridToWorld(cellId1);
    Vec2r wCoord2 = gridToWorld(cellId2);
    return (wCoord1-wCoord2).length();
}

std::array<Vec2r, 4> Grid2dUtility::corners(int pixelId) const
{
    Vec2r offset = gridToWorld(pixelId);
    std::array<Vec2r, 4> corners;
    corners[0] = offset;
    corners[1] = offset+Vec2r(h,0);
    corners[2] = offset+Vec2r(h,h);
    corners[3] = offset+Vec2r(0,h);
    return corners;
}

}//namespace mosaic

#endif
