#ifndef HOKUSAI_PARTICLE_CONTAINER_HPP
#define HOKUSAI_PARTICLE_CONTAINER_HPP

#include <vector>
#include <memory>

template<class ParticleT>
class ParticleContainer
{
public:
    ParticleContainer();
    ParticleContainer(const ParticleContainer<ParticleT>& other);
    ~ParticleContainer();

    std::size_t size() const;

    ParticleContainer& operator=(const ParticleContainer& other);

    ParticleT& operator[](std::size_t idx);
    const ParticleT& operator[](std::size_t idx) const;

    void add(const ParticleT& particle);

private:
    std::vector<ParticleT> m_particles;
};

template<class ParticleT>
class ParticleContainerPtr
{
public:
    ParticleContainerPtr();
    ParticleContainerPtr(ParticleContainerPtr<ParticleT>&& particles);
    ParticleContainerPtr(const ParticleContainerPtr<ParticleT>& particles);
    ~ParticleContainerPtr();

    std::size_t size() const;

    ParticleT& operator[](std::size_t idx);
    const ParticleT& operator[](std::size_t idx) const;

    ParticleContainerPtr& operator=(ParticleContainerPtr&& other);
    ParticleContainerPtr& operator=(ParticleContainerPtr other);

    void add(const ParticleT& particle);

    ParticleContainer<ParticleT>* get();
    const ParticleContainer<ParticleT>* get() const;

    void swap(ParticleContainerPtr& other);

private:
    std::shared_ptr< ParticleContainer<ParticleT> > m_particles;
};

#endif //HOKUSAI_PARTICLE_CONTAINER_HPP
