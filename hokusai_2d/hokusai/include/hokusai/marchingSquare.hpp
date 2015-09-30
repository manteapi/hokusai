#ifndef MARCHING_SQUARE_HPP
#define MARCHING_SQUARE_HPP

#include "utility.hpp"
#include "gridUtility.hpp"
#include <vector>
#include <array>

namespace hokusai
{

class Edge
{
public :
    Edge(){}
    Edge( Vec2d& _p1, Vec2d& _p2){ p1 = _p1; p2 = _p2;}
    ~Edge(){}
    Vec2d p1;
    Vec2d p2;
};

int sign(double a);

Vec2d vertexInterpolation(Vec2d& p1, Vec2d& p2, double value1, double value2, double isovalue);

//1------2
//|		 |
//|		 |
//4------3
void marchingSquare(std::vector<Edge>& edges, const std::array<double,4>& cellValue, const std::array<int,4>& cellSign,
                    const std::array<Vec2d,4>& cellVertex, const double isovalue);

void polygonize(std::vector<Edge>& edges, std::vector<double>& gridValue, const Grid2dUtility& gridInfo, const double isovalue);

} // namespace mosaic

#endif //MARCHING_SQUARE_HPP
