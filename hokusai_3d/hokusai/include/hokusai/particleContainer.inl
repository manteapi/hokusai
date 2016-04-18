#ifndef HOKUSAI_PARTICLE_CONTAINER_INL
#define HOKUSAI_PARTICLE_CONTAINER_INL

#include "particleContainer.hpp"

template<class ParticleT>
ParticleContainer<ParticleT>::ParticleContainer()
{

}

template<class ParticleT>
ParticleContainer<ParticleT>::ParticleContainer(const ParticleContainer<ParticleT>& other)
{
    m_particles = other.m_particles;
}

template<class ParticleT>
ParticleContainer<ParticleT>::~ParticleContainer()
{

}

template<class ParticleT>
std::size_t ParticleContainer<ParticleT>::size() const
{
    return m_particles.size();
}

template<class ParticleT>
void ParticleContainer<ParticleT>::add(const ParticleT& particle)
{
    m_particles.push_back(particle);
}

template<class ParticleT>
ParticleT& ParticleContainer<ParticleT>::operator[](std::size_t idx)
{
    return m_particles[idx];
};

template<class ParticleT>
const ParticleT& ParticleContainer<ParticleT>::operator[](std::size_t idx) const
{
    return m_particles[idx];
};

template<class ParticleT>
ParticleContainer<ParticleT>& ParticleContainer<ParticleT>::operator=(const ParticleContainer& other)
{
    m_particles = other.m_particles;
    return *this;
}

//===============================================

template<class ParticleT>
ParticleContainerPtr<ParticleT>::ParticleContainerPtr()
{
    m_particles = std::make_shared< ParticleContainer<ParticleT> >();
}

template<class ParticleT>
ParticleContainerPtr<ParticleT>::ParticleContainerPtr(ParticleContainerPtr<ParticleT>&& particles) :
    m_particles( std::move(particles.m_particles) )
{
    particles.m_particles = nullptr;
}

template<class ParticleT>
ParticleContainerPtr<ParticleT>::ParticleContainerPtr(const ParticleContainerPtr<ParticleT>& particles)
{
    m_particles = std::make_shared< ParticleContainer<ParticleT> >( *(particles.get()) ) ;
}

template<class ParticleT>
ParticleContainerPtr<ParticleT>::~ParticleContainerPtr()
{

}

template<class ParticleT>
std::size_t ParticleContainerPtr<ParticleT>::size() const
{
    return m_particles->size();
}

template<class ParticleT>
void ParticleContainerPtr<ParticleT>::add(const ParticleT& particle)
{
    m_particles->add(particle);
}

template<class ParticleT>
ParticleContainer<ParticleT>* ParticleContainerPtr<ParticleT>::get()
{
    return m_particles.get();
}

template<class ParticleT>
const ParticleContainer<ParticleT>* ParticleContainerPtr<ParticleT>::get() const
{
    return m_particles.get();
}

template<class ParticleT>
void ParticleContainerPtr<ParticleT>::swap(ParticleContainerPtr& other)
{
    std::swap(this->m_particles, other.m_particles);
}

template<class ParticleT>
ParticleT& ParticleContainerPtr<ParticleT>::operator[](std::size_t idx)
{
    return m_particles->operator[](idx);
};

template<class ParticleT>
const ParticleT& ParticleContainerPtr<ParticleT>::operator[](std::size_t idx) const
{
    return m_particles->operator[](idx);
};

template<class ParticleT>
ParticleContainerPtr<ParticleT>& ParticleContainerPtr<ParticleT>::operator=(ParticleContainerPtr&& other)
{
    m_particles = std::move(other.m_particles);
    return *this;
}

template<class ParticleT>
ParticleContainerPtr<ParticleT>& ParticleContainerPtr<ParticleT>::operator=(ParticleContainerPtr other)
{
    this->swap(other);
    return *this;
}

#endif //HOKUSAI_PARTICLE_CONTAINER_INL
