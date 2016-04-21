#ifndef HOKUSAI_SPATIAL_INDEX_INL
#define HOKUSAI_SPATIAL_INDEX_INL

#include "spatialIndex.hpp"

namespace hokusai
{

template<typename SpatialStructure>
SpatialIndex<SpatialStructure>::SpatialIndex(const StructureParams &params) :
    m_spatialStructure(params)
{}

template<typename SpatialStructure>
SpatialIndex<SpatialStructure>::~SpatialIndex()
{

}

template<typename SpatialStructure>
void SpatialIndex<SpatialStructure>::buildIndex( const std::vector<Vec3r>& data )
{
    m_spatialStructure.buildIndex(data);
}

template<typename SpatialStructure>
void SpatialIndex<SpatialStructure>::radiusSearch( std::vector<int>& indices, std::vector<HReal>& dists, const Vec3r& query, const HReal& radius)
{
    m_spatialStructure.radiusSearch(indices, dists, query, radius);
}

template<typename SpatialStructure>
void SpatialIndex<SpatialStructure>::radiusSearch(std::vector< std::vector<int> >& indices, std::vector< std::vector<HReal> > &dists, const std::vector<Vec3r>& query, const HReal& radius)
{
    for(size_t i=0; i<query.size(); ++i)
    {
        m_spatialStructure.radiusSearch(indices[i], dists[i], query[i], radius);
    }
}

} //namespace hokusai

#endif //HOKUSAI_SPATIAL_INDEX_INL
