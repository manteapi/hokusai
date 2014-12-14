#ifndef GRID_UTILITY_CPP
#define GRID_UTILITY_CPP

#include "../include/hokusai/gridUtility.hpp"

GridUtility::GridUtility()
{
    offset = Vec3f(0,0,0);
    scale = Vec3f(0,0,0);
    dimension = Vec3i(0,0,0);
    h = 0.0;
}

int GridUtility::width(){ return dimension[0]; }
int GridUtility::height(){ return dimension[1]; }
int GridUtility::depth(){ return dimension[2]; }
float GridUtility::spacing(){ return h; }

GridUtility::GridUtility(const Vec3f& _offset, const Vec3f& _scale, const float& _spacing)
{
    offset = _offset;
    scale = _scale;
    h = _spacing;
    for(int i=0; i<3; ++i)
    {
        dimension[i] = std::floor(scale[i]/h);
    }
}

GridUtility::GridUtility(const Vec3f& _offset, const Vec3i& _dimension, const float& _spacing)
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

Vec3f GridUtility::gridToWorld(const Vec3i& v)
{
    Vec3f result;
    for(int i=0; i<3; ++i)
        result[i] = offset[i] + h*v[i];
    return result;
}
Vec3i GridUtility::worldToGrid(const Vec3f& v)
{
    Vec3i result;
    for(int i=0; i<3; ++i)
        result[i] = std::floor((v[i]-offset[i])/h);
    return result;
}
int GridUtility::cellId(const Vec3f& v)
{
    Vec3i gridCoord = worldToGrid(v);
    return gridCoord[0] + gridCoord[1]*std::floor(scale[0]/h) + gridCoord[2]*std::floor(scale[0]/h)*std::floor(scale[1]/h);
}

int GridUtility::cellId(const Vec3i& v)
{
    return v[0] + v[1]*std::floor(scale[0]/h) + v[2]*std::floor(scale[0]/h)*std::floor(scale[1]/h);
}

int GridUtility::cellId(int i, int j, int k)
{
    return i + j*std::floor(scale[0]/h) + k*std::floor(scale[0]/h)*std::floor(scale[1]/h);
}

#endif
