#ifndef HOKUSAI_PARTICLE_CONTAINER_HPP
#define HOKUSAI_PARTICLE_CONTAINER_HPP

#include <vector>

template<class ParticleT>
class ParticleContainer
{
public:
    ParticleContainer();
    ~ParticleContainer();

    std::size_t size() const;

    ParticleContainer& operator=(const ParticleContainer& other);

    ParticleT& operator[](std::size_t idx);
    const ParticleT& operator[](std::size_t idx) const;

    void add(const ParticleT& particle);

private:
    std::vector<ParticleT> m_particles;
};

#endif //HOKUSAI_PARTICLE_CONTAINER_HPP
