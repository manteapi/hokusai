#ifndef IMPLICIT_FUNCTION_CPP
#define IMPLICIT_FUNCTION_CPP

#include "./../include/hokusai/implicitFunction.hpp"

namespace hokusai
{

void computeScalarField(std::vector<double>& scalarField, Grid2dUtility& gridInfo, System& fluid, const double& resolution, const double& initialValue)
{
    //Compute scalar field dimensions
    Vec2d minBB(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    Vec2d maxBB(-std::numeric_limits<double>::max(),-std::numeric_limits<double>::max());
    for(size_t i=0; i<fluid.particles.size(); ++i)
    {
        const Vec2d& px = fluid.particles[i].x;
        minBB[0] = (px[0] < minBB[0]) ? px[0] : minBB[0];
        minBB[1] = (px[1] < minBB[1]) ? px[1] : minBB[1];
        maxBB[0] = (px[0] > maxBB[0]) ? px[0] : maxBB[0];
        maxBB[1] = (px[1] > maxBB[1]) ? px[1] : maxBB[1];
    }
    double radius = fluid.getSmoothingRadius();
    Vec2d offset = minBB - Vec2d(2.0*radius, 2.0*radius);
    Vec2d scale = (maxBB-minBB) + Vec2d(4.0*radius, 4.0*radius);
    if(fluid.particles.size()>0)
        gridInfo = Grid2dUtility(offset, scale, resolution);
    else
        gridInfo = Grid2dUtility();

    //Compute scalar field particle neighbors
    double implicitFunctionRadius = radius;
    std::vector< std::vector<int> > gridParticleNeighbors;
    gridParticleNeighbors.resize(gridInfo.size());
#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<gridInfo.size(); ++i)
    {
        Vec2d gridWCoord = gridInfo.gridToWorld(i);
        fluid.getNearestFluidNeighbor(gridParticleNeighbors[i], gridWCoord, implicitFunctionRadius);
    }

    //Compute scalar field values
    scalarField.resize(gridInfo.size());
    std::fill(scalarField.begin(), scalarField.end(), initialValue);
#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<gridInfo.size(); ++i)
    {
        Vec2d gridWCoord = gridInfo.gridToWorld(i);
        std::vector<int> neighbors = gridParticleNeighbors[i];
        if(neighbors.size()>0)
            scalarField[i] = 0.0;
        for(const int & j : neighbors)
        {
            Vec2d particleCoord = fluid.particles[j].x;
            scalarField[i] += implicit_sphere(gridWCoord, particleCoord, implicitFunctionRadius);
        }
    }
}

double implicit_sphere(const Vec2d & x1, const Vec2d & x2, const double radius)
{
    return ( (x1-x2).squaredNorm()-radius*radius );
}

} //namespace hokusai

#endif //IMPLICIT_FUNCTION_CPP
