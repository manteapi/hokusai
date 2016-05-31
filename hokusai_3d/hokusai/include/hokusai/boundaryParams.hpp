#ifndef HOKUSAI_BOUNDARY_PARAMS_HPP
#define HOKUSAI_BOUNDARY_PARAMS_HPP

#include "common.hpp"
#include "kernel.hpp"
#include <memory>

namespace hokusai
{

class BoundaryParams
{
public:
    BoundaryParams();
    BoundaryParams(HReal radius, HReal adhesion=0.001, HReal friction=1.0);
    ~BoundaryParams();

    HReal& radius();
    const HReal& radius() const;

    HReal& adhesion();
    const HReal& adhesion() const;

    HReal& friction();
    const HReal& friction() const;

    HReal& soundSpeed();
    const HReal& soundSpeed() const;

    BoundaryKernel& boundaryKernel();
    const BoundaryKernel& boundaryKernel() const;

private:
    HReal m_radius;
    HReal m_adhesion;
    HReal m_friction;
    HReal m_soundSpeed;
    BoundaryKernel m_bKernel;
};

typedef std::shared_ptr<BoundaryParams> BoundaryParamsPtr;

}//namespace hokusai

#endif//HOKUSAI_BOUNDARY_PARAMS_HPP
