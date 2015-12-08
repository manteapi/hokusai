#ifndef GRID_UTILITY_CPP
#define GRID_UTILITY_CPP

#include "../include/hokusai/gridUtility.hpp"

namespace hokusai
{

Grid2dUtility::Grid2dUtility()
{
    offset = Vec2d(0,0);
    scale = Vec2d(0,0);
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

Grid2dUtility::Grid2dUtility(const Vec2d& _offset, const Vec2d& _scale, const double& _spacing)
{
    offset = _offset;
    scale = _scale;
    h = _spacing;
    for(int i=0; i<scale.size(); ++i)
    {
        dimension[i] = std::floor(scale[i]/h);
    }
}

Grid2dUtility::Grid2dUtility(const Vec2d& _offset, const Vec2i& _dimension, const double& _spacing)
{
    offset = _offset;
    dimension = _dimension;
    h = _spacing;
    for(int i=0; i<scale.size(); ++i)
    {
        scale[i] = dimension[i]*h;
    }
}

void Grid2dUtility::update(const Vec2d& _offset, const Vec2i& _dimension, const double& _spacing)
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

void Grid2dUtility::update(const Vec2d& _offset, const Vec2d& _scale, const double& _spacing)
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

void Grid2dUtility::init(const Vec2d& _offset, const Vec2d &_scale, const double& _spacing)
{
    offset = _offset;
    scale = _scale;
    h = _spacing;
    for(int i=0; i<2; ++i)
    {
        dimension[i] = std::floor(scale[i]/h);
    }
}

void Grid2dUtility::init(const Vec2d& _offset, const Vec2i& _dimension, const double& _spacing)
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

void Grid2dUtility::info()
{
    std::cout << "Dimension: " << width() << "x" << height() << " = " << size() << std::endl;
    std::cout << "Spacing : " << spacing() << std::endl;
    std::cout << "Offset : " << offset[0] << ", " << offset[1] << std::endl;
    std::cout << "Scale : " << scale[0] << ", " << scale[1] << std::endl;
    std::cout << "Max : " << offset[0]+scale[0] << ", " << offset[1]+scale[1] << std::endl;
}

}//namespace mosaic

#endif
