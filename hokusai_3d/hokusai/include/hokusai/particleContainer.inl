#ifndef HOKUSAI_PARTICLE_CONTAINER_INL
#define HOKUSAI_PARTICLE_CONTAINER_INL

#include "particleContainer.hpp"

template<class ParticleT>
ParticleContainer<ParticleT>::ParticleContainer()
{

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

#endif //HOKUSAI_PARTICLE_CONTAINER_INL
