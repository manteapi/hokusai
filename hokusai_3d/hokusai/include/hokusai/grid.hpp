#ifndef HOKUSAI_GRID_HPP
#define HOKUSAI_GRID_HPP

#include "common.hpp"
#include "gridUtility.hpp"

namespace hokusai
{

class GridParams
{
public:
    GridParams(const HReal& spacing);
    ~GridParams();
    const HReal& spacing() const;
private:
    HReal m_spacing;
};

class Grid
{
public:
    typedef GridParams Params;
    Grid( const GridParams& params);
    ~Grid();
    void buildIndex( const std::vector<Vec3r>& data );
    void radiusSearch( std::vector<int>& indices, std::vector<HReal>& dists, const Vec3r& query, const HReal& radius);
    GridUtility m_gridInfo;
    std::vector< std::vector<int> > m_gridData;
};

}//namespace hokusai

#endif
