#include "./../include/hokusai/boundaryParams.hpp"

namespace hokusai
{

BoundaryParams::BoundaryParams()
{

}

BoundaryParams::BoundaryParams(const HReal& smoothingRadius, const HReal& adhesion, const HReal& friction)
{
    m_smoothingRadius = smoothingRadius;
    m_adhesion = adhesion;
    m_friction = friction;
    HReal eta = 0.01;
    HReal H = 0.1;
    HReal vf = sqrt( 2*9.81*H );
    m_soundSpeed = vf/(sqrt(eta));
    m_boundaryKernel = BoundaryKernel(m_smoothingRadius, m_soundSpeed);
}

BoundaryParams::BoundaryParams(const BoundaryParams& params)
{
    m_smoothingRadius = params.smoothingRadius();
    m_adhesion = params.adhesion();
    m_friction = params.friction();
    m_boundaryKernel = params.boundaryKernel();
}

BoundaryParams::~BoundaryParams()
{

}

HReal& BoundaryParams::adhesion()
{
    return m_adhesion;
}

const HReal& BoundaryParams::adhesion() const
{
    return m_adhesion;
}

HReal& BoundaryParams::friction()
{
    return m_friction;
}

const HReal& BoundaryParams::friction() const
{
    return m_friction;
}

HReal& BoundaryParams::smoothingRadius()
{
    return m_smoothingRadius;
}

const HReal& BoundaryParams::smoothingRadius() const
{
    return m_smoothingRadius;
}

HReal& BoundaryParams::soundSpeed()
{
    return m_soundSpeed;
}

const HReal& BoundaryParams::soundSpeed() const
{
    return m_soundSpeed;
}

BoundaryKernel& BoundaryParams::boundaryKernel()
{
    return m_boundaryKernel;
}

const BoundaryKernel& BoundaryParams::boundaryKernel() const
{
    return m_boundaryKernel;
}

} //namespace hokusai
