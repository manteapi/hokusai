#ifndef HOKUSAI_SPATIAL_INDEX_HPP
#define HOKUSAI_SPATIAL_INDEX_HPP

#include "common.hpp"
#include <memory>

namespace hokusai
{

template<typename SpatialStructure>
class SpatialIndex
{
    typedef typename SpatialStructure::Params StructureParams;

public:
    SpatialIndex(const StructureParams& params);
    ~SpatialIndex();
    SpatialStructure m_spatialStructure;
    void buildIndex( const std::vector<Vec3r>& data );
    void radiusSearch( std::vector<int>& indices, std::vector<HReal>& dists, const Vec3r& query, const HReal& radius);
    void radiusSearch( std::vector< std::vector<int> >& indices, std::vector< std::vector<HReal> >& dists, const std::vector<Vec3r>& query, const HReal& radius);
};

} //namespace hokusai

#endif //HOKUSAI_SPATIAL_INDEX_HPP
