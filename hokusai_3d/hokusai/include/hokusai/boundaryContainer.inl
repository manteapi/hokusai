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

//===============================================

template<class BoundaryT>
BoundaryContainerPtr<BoundaryT>::BoundaryContainerPtr()
{
    m_boundaries = std::make_shared< BoundaryContainer<BoundaryT> >();
}

template<class BoundaryT>
BoundaryContainerPtr<BoundaryT>::BoundaryContainerPtr(BoundaryContainerPtr<BoundaryT>&& boundaries) :
    m_boundaries( std::move(boundaries.m_boundaries) )
{
    boundaries.m_boundaries = nullptr;
}

template<class BoundaryT>
BoundaryContainerPtr<BoundaryT>::BoundaryContainerPtr(const BoundaryContainerPtr<BoundaryT>& boundaries)
{
    m_boundaries = std::make_shared< BoundaryContainer<BoundaryT> >( *(boundaries.get()) ) ;
}

template<class BoundaryT>
BoundaryContainerPtr<BoundaryT>::~BoundaryContainerPtr()
{

}

template<class BoundaryT>
std::size_t BoundaryContainerPtr<BoundaryT>::size() const
{
    return m_boundaries->size();
}

template<class BoundaryT>
void BoundaryContainerPtr<BoundaryT>::add(const BoundaryT& particle)
{
    m_boundaries->add(particle);
}

template<class BoundaryT>
BoundaryContainer<BoundaryT>* BoundaryContainerPtr<BoundaryT>::get()
{
    return m_boundaries.get();
}

template<class BoundaryT>
const BoundaryContainer<BoundaryT>* BoundaryContainerPtr<BoundaryT>::get() const
{
    return m_boundaries.get();
}

template<class BoundaryT>
void BoundaryContainerPtr<BoundaryT>::swap(BoundaryContainerPtr& other)
{
    std::swap(this->m_boundaries, other.m_boundaries);
}

template<class BoundaryT>
BoundaryT& BoundaryContainerPtr<BoundaryT>::operator[](std::size_t idx)
{
    return m_boundaries->operator[](idx);
};

template<class BoundaryT>
const BoundaryT& BoundaryContainerPtr<BoundaryT>::operator[](std::size_t idx) const
{
    return m_boundaries->operator[](idx);
};

template<class BoundaryT>
BoundaryContainerPtr<BoundaryT>& BoundaryContainerPtr<BoundaryT>::operator=(BoundaryContainerPtr&& other)
{
    m_boundaries = std::move(other.m_boundaries);
    return *this;
}

template<class BoundaryT>
BoundaryContainerPtr<BoundaryT>& BoundaryContainerPtr<BoundaryT>::operator=(BoundaryContainerPtr other)
{
    this->swap(other);
    return *this;
}

#endif //HOKUSAI_BOUNDARY_CONTAINER_INL
