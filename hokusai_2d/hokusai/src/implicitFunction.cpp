#ifndef IMPLICIT_FUNCTION_CPP
#define IMPLICIT_FUNCTION_CPP

#include "./../include/hokusai/implicitFunction.hpp"

namespace hokusai
{

void computeScalarField(std::vector<double>& scalarField, Grid2dUtility& gridInfo, System& fluid,
                        const double& resolution, const double& initialValue, const double& radius)
{
    //Compute scalar field dimensions
    Vec2d minBB(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    Vec2d maxBB(-std::numeric_limits<double>::max(),-std::numeric_limits<double>::max());
    for(size_t i=0; i<fluid.particles.size(); ++i)
    {
        const Vec2d& px = fluid.particles[i].x;
        for(size_t j=0; j<2; ++j)
        {
            minBB[j] = (px[j] < minBB[j]) ? px[j] : minBB[j];
            maxBB[j] = (px[j] > maxBB[j]) ? px[j] : maxBB[j];
        }
    }
    if(resolution>radius)
    {
        std::cout << "Grid resolution should be at least inferior to sph radius" << std::endl;
        return;
    }
    Vec2d scale = (maxBB-minBB) + Vec2d(16.0*std::max(radius,resolution), 16.0*std::max(radius,resolution));
    Vec2d offset = minBB - Vec2d(8.0*std::max(radius,resolution), 8.0*std::max(radius,resolution));
    if(fluid.particles.size()>0)
        gridInfo = Grid2dUtility(offset, scale, resolution);
    else
        gridInfo = Grid2dUtility();

    scalarField.resize(gridInfo.size());
    std::fill(scalarField.begin(), scalarField.end(), initialValue);
    for(size_t p=0; p<fluid.particles.size(); ++p)
    {
        Vec2d & wp = fluid.particles[p].x;
        Vec2d wMinBB = wp - Vec2d(radius, radius);
        Vec2d wMaxBB = wp + Vec2d(radius, radius);
        Vec2i gMinBB = gridInfo.worldToGrid(wMinBB);
        Vec2i gMaxBB = gridInfo.worldToGrid(wMaxBB);
        for(int i=gMinBB[0]; i<=gMaxBB[0]; ++i)
        {
            for(int j=gMinBB[1]; j<=gMaxBB[1]; ++j)
            {
                int cellId = gridInfo.cellId(i,j);
                if(gridInfo.isInside(cellId))
                {
                    Vec2d wCoord = gridInfo.gridToWorld(cellId);
                    ///Sphere implicit function
                    //double implicitFunctionValue = implicit_sphere(wCoord, fluid.particles[p].x, worldRadius);
                    //scalarField[cellId] = std::min( std::min(0.0, implicitFunctionValue), scalarField[cellId]);
                    ///Metaball implicit function
                    double implicitFunctionValue = implicit_metaball(wCoord, fluid.particles[p].x,radius);
                    scalarField[cellId] += implicitFunctionValue;
                }
            }
        }
    }
}

double implicit_sphere(const Vec2d & x1, const Vec2d & x2, const double radius)
{
    return ( (x1-x2).squaredNorm()-radius*radius );
}

double implicit_metaball(const Vec2d & x1, const Vec2d & x2, const double radius)
{
    double distance=(x1-x2).norm();
    return distance<radius ? 1.0/(1.0+distance*distance) : 0.0;
}

} //namespace hokusai

#endif //IMPLICIT_FUNCTION_CPP
