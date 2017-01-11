#ifndef HOKUSAI_BOUNDARY_PARAMS_HPP
#define HOKUSAI_BOUNDARY_PARAMS_HPP

#include "common.hpp"
#include "kernel.hpp"

namespace hokusai
{

class BoundaryParams
{
public:
    BoundaryParams();
    BoundaryParams(const HReal& smoothingRadius, const HReal& adhesion, const HReal& friction);
    BoundaryParams(const BoundaryParams& params);
    ~BoundaryParams();
    HReal& adhesion();
    const HReal& adhesion() const;
    HReal& friction();
    const HReal& friction() const;
    HReal& smoothingRadius();
    const HReal& smoothingRadius() const;
    BoundaryKernel& boundaryKernel();
    const BoundaryKernel& boundaryKernel() const;
    HReal& soundSpeed();
    const HReal& soundSpeed() const;
private:
    HReal m_soundSpeed;
    HReal m_adhesion;
    HReal m_friction;
    HReal m_smoothingRadius;
    BoundaryKernel m_boundaryKernel;
};

} //namespace hokusai

#endif //HOKUSAI_BOUNDARY_PARAMS_HPP
