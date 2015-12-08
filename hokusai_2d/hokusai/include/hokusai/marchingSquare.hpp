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
    Edge(const Edge& e){ p1 = e.p1; p2 = e.p2;}
    ~Edge(){}
    Vec2d p1;
    Vec2d p2;
};

//bool operator== (const Edge& e0, const Edge& e1)
//{
//    if((e0.p1 - e1.p1).norm() < 1e-10 && (e0.p2 - e1.p2).norm() < 1e-10)
//    {
//        return true;
//    }
//    return false;
//}

inline int sign(double a)
{
    if( a < 0){return -1;}
    else {return 1;}
}

inline Vec2d vertexInterpolation(const Vec2d& p1, const Vec2d& p2, double value1, double value2, double isovalue)
{
    Vec2d interp(0.0,0.0);
    double weight;
    ///If value are too close to each other then I used to snap it to the closest point.
    ///However in practice it might not snap the point to the expected point.
    ///For now I unactivated the test (21/10/2015)
    double epsilon = 0;//std::numeric_limits<double>::epsilon();
    if( abs(isovalue-value1) < epsilon ){ interp = p1; return interp;}
    if( abs(isovalue-value2) < epsilon ){ interp = p2; return interp;}
    weight = (isovalue-value1)/(value2-value1);
    interp = p1 + weight * ( p2 - p1 );
    return interp;
}

//1------2
//|		 |
//|		 |
//4------3
//For a cell, compute edges
void marchingSquare(std::vector<Edge>& edges, const std::array<double,4>& cellValue, const std::array<int,4>& cellSign,
                    const std::array<Vec2d,4>& cellVertex, const double &isovalue);

//For each cell of a grid, compute edges
void polygonize(std::vector<Edge>& edges, std::vector<double>& gridValue, const Grid2dUtility& gridInfo, const double& isovalue);

} // namespace mosaic

#endif //MARCHING_SQUARE_HPP
