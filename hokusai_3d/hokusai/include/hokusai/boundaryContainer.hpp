#ifndef HOKUSAI_BOUNDARY_CONTAINER_HPP
#define HOKUSAI_BOUNDARY_CONTAINER_HPP

#include <vector>

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

#endif //HOKUSAI_BOUNDARY_CONTAINER_HPP
