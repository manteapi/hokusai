#ifndef MARCHING_SQUARE_CPP
#define MARCHING_SQUARE_CPP

#include "./../include/hokusai/marchingSquare.hpp"

namespace hokusai
{

int sign(double a)
{
    if( a < 0){return -1;}
    else {return 1;}
}

Vec2d vertexInterpolation(const Vec2d& p1, const Vec2d& p2, double value1, double value2, double isovalue)
{
    Vec2d interp;
    double weight;
    double epsilon = 1e-4;
    if( abs(isovalue-value1) < epsilon ){ interp = p1; return interp;}
    if( abs(isovalue-value2) < epsilon ){ interp = p2; return interp;}
    weight = (isovalue-value1)/(value2-value1);
    interp = p1 + weight * ( p2 - p1 );
    return interp;
}

void marchingSquare(std::vector<Edge>& edges, const std::array<double,4>& cellValue, const std::array<int,4>& cellSign,
                    const std::array<Vec2d,4>& cellVertex, const double isovalue)
{
    //Case 1 : all same sign = no edges
    if( (cellSign[0] == cellSign[1]) && (cellSign[1] == cellSign[2]) && (cellSign[2]==cellSign[3]) ){/*nothing to do*/ return;}

    //Case 2 : one sign only (8 possibilities)
    //1+-
    if(	( (cellSign[0] == -1) && (cellSign[1] == 1) && (cellSign[2] == 1) && (cellSign[3] == 1) ) ||
            ( (cellSign[0] == 1) && (cellSign[1] == -1) && (cellSign[2] == -1) && (cellSign[3] == -1) ) )
    {
        Vec2d p1 = vertexInterpolation(cellVertex[3], cellVertex[0], cellValue[3], cellValue[0], isovalue);
        Vec2d p2 = vertexInterpolation(cellVertex[0], cellVertex[1], cellValue[0], cellValue[1], isovalue);
        edges.push_back(Edge(p1,p2));
        return;
    }
    //2+-
    if(	( (cellSign[0] == 1) && (cellSign[1] == -1) && (cellSign[2] == 1) && (cellSign[3] == 1) ) ||
            ( (cellSign[0] == -1) && (cellSign[1] == 1) && (cellSign[2] == -1) && (cellSign[3] == -1) ) )
    {
        Vec2d p1 = vertexInterpolation(cellVertex[0], cellVertex[1], cellValue[0], cellValue[1], isovalue);
        Vec2d p2 = vertexInterpolation(cellVertex[1], cellVertex[2], cellValue[1], cellValue[2], isovalue);
        edges.push_back(Edge(p1,p2));
        return;
    }
    //3+-
    if(	( (cellSign[0] == 1) && (cellSign[1] == 1) && (cellSign[2] == -1) && (cellSign[3] == 1) ) ||
            ( (cellSign[0] == -1) && (cellSign[1] == -1) && (cellSign[2] == 1) && (cellSign[3] == -1) ) )
    {
        Vec2d p1 = vertexInterpolation(cellVertex[1], cellVertex[2], cellValue[1], cellValue[2], isovalue);
        Vec2d p2 = vertexInterpolation(cellVertex[2], cellVertex[3], cellValue[2], cellValue[3], isovalue);
        edges.push_back(Edge(p1,p2));
        return;
    }
    //4+-
    if(	( (cellSign[0] == 1) && (cellSign[1] == 1) && (cellSign[2] == 1) && (cellSign[3] == -1) ) ||
            ( (cellSign[0] == -1) && (cellSign[1] == -1) && (cellSign[2] == -1) && (cellSign[3] == 1) ) )
    {
        Vec2d p1 = vertexInterpolation(cellVertex[2], cellVertex[3], cellValue[2], cellValue[3], isovalue);
        Vec2d p2 = vertexInterpolation(cellVertex[3], cellVertex[0], cellValue[3], cellValue[0], isovalue);
        edges.push_back(Edge(p1,p2));
        return;
    }

    //Case 3 : two adjacent same sign ( 4 possibilities)
    //horizontal
    if( ( (cellSign[0]==-1) && (cellSign[1]==-1) && (cellSign[2]==1) && (cellSign[3]==1) ) ||
            ( (cellSign[0]==1) && (cellSign[1]==1) && (cellSign[2]==-1) && (cellSign[3]==-1) ) )
    {
        Vec2d p1 = vertexInterpolation(cellVertex[3], cellVertex[0], cellValue[3], cellValue[0], isovalue);
        Vec2d p2 = vertexInterpolation(cellVertex[1], cellVertex[2], cellValue[1], cellValue[2], isovalue);
        edges.push_back(Edge(p1,p2));
        return;
    }
    //vertical
    if( ( (cellSign[0]==-1) && (cellSign[1]==1) && (cellSign[2]==1) && (cellSign[3]==-1) ) ||
            ( (cellSign[0]==1) && (cellSign[1]==-1) && (cellSign[2]==-1) && (cellSign[3]==1) ) )
    {
        Vec2d p1 = vertexInterpolation(cellVertex[0], cellVertex[1], cellValue[0], cellValue[1], isovalue);
        Vec2d p2 = vertexInterpolation(cellVertex[2], cellVertex[3], cellValue[2], cellValue[3], isovalue);
        edges.push_back(Edge(p1,p2));
        return;
    }

    //Case 4 : two diagonal same sign (4 possibilite)
    // - + - +
    //Compute center value
    double centerValue = 0.25*(cellValue[0] + cellValue[1] + cellValue[2] + cellValue[3]);
    int sign5 = sign(centerValue-isovalue);

    if( ( (cellSign[0]==-1) && (cellSign[1]==1) && (cellSign[2]==-1) && (cellSign[3]==1) && (sign5==1) ) ||
            ( (cellSign[0]==1) && (cellSign[1]==-1) && (cellSign[2]==1) && (cellSign[3]==-1) && (sign5==-1) ) )
    {
        Vec2d p1 = vertexInterpolation(cellVertex[3], cellVertex[0], cellValue[3], cellValue[0], isovalue);
        Vec2d p2 = vertexInterpolation(cellVertex[0], cellVertex[1], cellValue[0], cellValue[1], isovalue);
        edges.push_back(Edge(p1,p2));

        Vec2d p3 = vertexInterpolation(cellVertex[1], cellVertex[2], cellValue[1], cellValue[2], isovalue);
        Vec2d p4 = vertexInterpolation(cellVertex[2], cellVertex[3], cellValue[2], cellValue[3], isovalue);
        edges.push_back(Edge(p3,p4));
        return;
    }
    if( ( (cellSign[0]==-1) && (cellSign[1]==1) && (cellSign[2]==-1) && (cellSign[3]==1) && (sign5==-1) ) ||
            ( (cellSign[0]==1) && (cellSign[1]==-1) && (cellSign[2]==1) && (cellSign[3]==-1) && (sign5==1) ) )
    {
        Vec2d p1 = vertexInterpolation(cellVertex[0], cellVertex[1], cellValue[0], cellValue[1], isovalue);
        Vec2d p2 = vertexInterpolation(cellVertex[1], cellVertex[2], cellValue[1], cellValue[2], isovalue);
        edges.push_back(Edge(p1,p2));

        Vec2d p3 = vertexInterpolation(cellVertex[2], cellVertex[3], cellValue[2], cellValue[3], isovalue);
        Vec2d p4 = vertexInterpolation(cellVertex[3], cellVertex[0], cellValue[3], cellValue[0], isovalue);
        edges.push_back(Edge(p3,p4));
        return;
    }
}


void polygonize(std::vector<Edge>& edges, std::vector<double>& gridValue, const Grid2dUtility& gridInfo, const double isovalue)
{
    edges.clear();
    for(int i=0; i<gridInfo.width()-1; ++i)
    {
        for(int j=0; j<gridInfo.height()-1; ++j)
        {
            std::array<int,4> cellId = {gridInfo.cellId(i,j+1), gridInfo.cellId(i+1,j+1), gridInfo.cellId(i+1,j), gridInfo.cellId(i,j)};
            std::array<double,4> cellValue = {gridValue[cellId[0]], gridValue[cellId[1]], gridValue[cellId[2]], gridValue[cellId[3]]};
            std::array<int,4> cellSign = {sign(cellValue[0] - isovalue), sign(cellValue[1] - isovalue), sign(cellValue[2] - isovalue), sign(cellValue[3] - isovalue)};
            std::array<Vec2d,4> cellVertex = {gridInfo.gridToWorld(cellId[0]), gridInfo.gridToWorld(cellId[1]), gridInfo.gridToWorld(cellId[2]), gridInfo.gridToWorld(cellId[3])};
            marchingSquare(edges, cellValue, cellSign, cellVertex, isovalue);
        }
    }
}


} //namespace hokusai

#endif //MARCHING_SQUARE_CPP
