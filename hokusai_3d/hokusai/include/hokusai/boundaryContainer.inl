#ifndef HOKUSAI_BOUNDARY_CONTAINER_INL
#define HOKUSAI_BOUNDARY_CONTAINER_INL

#include "boundaryContainer.hpp"

template<class BoundaryT>
BoundaryContainer<BoundaryT>::BoundaryContainer()
{

}

template<class BoundaryT>
BoundaryContainer<BoundaryT>::~BoundaryContainer()
{

}

template<class BoundaryT>
std::size_t BoundaryContainer<BoundaryT>::size() const
{
    return m_boundaries.size();
}

template<class BoundaryT>
void BoundaryContainer<BoundaryT>::add(const BoundaryT& particle)
{
    m_boundaries.push_back(particle);
}

template<class BoundaryT>
BoundaryT& BoundaryContainer<BoundaryT>::operator[](std::size_t idx)
{
    return m_boundaries[idx];
};

template<class BoundaryT>
const BoundaryT& BoundaryContainer<BoundaryT>::operator[](std::size_t idx) const
{
    return m_boundaries[idx];
};

template<class BoundaryT>
BoundaryContainer<BoundaryT>& BoundaryContainer<BoundaryT>::operator=(const BoundaryContainer& other)
{
    m_boundaries = other.m_boundaries;
    return *this;
}

#endif //HOKUSAI_BOUNDARY_CONTAINER_INL
