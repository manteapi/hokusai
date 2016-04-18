#ifndef HOKUSAI_BOUNDARY_CONTAINER_HPP
#define HOKUSAI_BOUNDARY_CONTAINER_HPP

#include <vector>
#include <memory>

template<class BoundaryT>
class BoundaryContainer
{
public:
    BoundaryContainer();
    ~BoundaryContainer();

    std::size_t size() const;

    BoundaryContainer& operator=(const BoundaryContainer& other);

    BoundaryT& operator[](std::size_t idx);
    const BoundaryT& operator[](std::size_t idx) const;

    void add(const BoundaryT& particle);

private:
    std::vector<BoundaryT> m_boundaries;
};

template<class BoundaryT>
class BoundaryContainerPtr
{
public:
    BoundaryContainerPtr();
    BoundaryContainerPtr(BoundaryContainerPtr<BoundaryT>&& boundaries);
    BoundaryContainerPtr(const BoundaryContainerPtr<BoundaryT>& boundaries);
    ~BoundaryContainerPtr();

    std::size_t size() const;

    BoundaryT& operator[](std::size_t idx);
    const BoundaryT& operator[](std::size_t idx) const;

    BoundaryContainerPtr& operator=(BoundaryContainerPtr&& other);
    BoundaryContainerPtr& operator=(BoundaryContainerPtr other);

    void add(const BoundaryT& Boundary);

    BoundaryContainer<BoundaryT>* get();
    const BoundaryContainer<BoundaryT>* get() const;

    void swap(BoundaryContainerPtr& other);

private:
    std::shared_ptr< BoundaryContainer<BoundaryT> > m_boundaries;
};

#endif //HOKUSAI_BOUNDARY_CONTAINER_HPP
