#include "./../include/hokusai/grid.hpp"

namespace hokusai
{

GridParams::GridParams(const HReal& spacing)
{
    m_spacing = spacing;
}

GridParams::~GridParams()
{

}

const HReal& GridParams::spacing() const
{
    return m_spacing;
}

Grid::Grid(const GridParams &params)
{
    Vec3r offset(0.0), scale(0.0);
    m_gridInfo = GridUtility(offset, scale, params.spacing());
}

Grid::~Grid()
{

}

void Grid::buildIndex( const std::vector<Vec3r>& data )
{
    Vec3r minBB(std::numeric_limits<HReal>::max()), maxBB(-std::numeric_limits<HReal>::max());
    for(size_t i=0; i<data.size(); ++i)
    {
        for(int j=0; j<3; ++j)
        {
            minBB[j] = std::min(data[i][j], minBB[j]);
            maxBB[j] = std::max(data[i][j], maxBB[j]);
        }
    }
    HReal spacing = m_gridInfo.spacing();
    Vec3r offset = minBB;
    Vec3r scale = maxBB-minBB;
    m_gridInfo = GridUtility( offset-Vec3r(spacing), scale+Vec3r(2.0*spacing), spacing );
}

void Grid::radiusSearch( std::vector<int>& indices, std::vector<HReal>& dists, const Vec3r& query, const HReal& radius)
{
/*
 m_gridInfo.get27Neighbors(indices, query, radius);

    for(size_t i=0; i<neighborCell.size(); ++i)
    {
        std::vector<int>& bNeighborCell = m_boundaryGrid[neighborCell[i]];
        for(size_t j=0; j<bNeighborCell.size(); ++j)
        {
            int bParticleId = m_boundaryGrid[neighborCell[i]][j];
            BoundaryIISPH& bParticle = m_boundaries[bParticleId];
            Vec3r d = bParticle.x-p.x;
            if( d.lengthSquared()<radius*radius )
                p.boundaryNeighbor.push_back(bParticleId);
        }

        std::vector<int>& fNeighborCell = m_fluidGrid[neighborCell[i]];
        for(size_t j=0; j< fNeighborCell.size(); ++j)
        {
            int fParticleId = m_fluidGrid[neighborCell[i]][j];
            ParticleIISPH& fParticle = m_particles[fParticleId];
            Vec3r d = fParticle.x-p.x;
            if( d.lengthSquared()<radius*radius)
                p.fluidNeighbor.push_back(fParticleId);
        }
    }
    */
}

}//namespace hokusai
