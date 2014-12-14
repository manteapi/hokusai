#ifndef GRID_UTILITY_HPP
#define GRID_UTILITY_HPP

#include "Vec.hpp"

typedef Vec3<float> Vec3f;
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

        Vec3f gridToWorld(const Vec3i& v);
        Vec3i worldToGrid(const Vec3f& v);
        int cellId(const Vec3f& v);
        int cellId(const Vec3i& v);
        int cellId(int i, int j, int k);
        int width();
        int height();
        int depth();
        float spacing();
};

#endif
