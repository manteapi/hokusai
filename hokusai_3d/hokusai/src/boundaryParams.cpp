#include "./../include/hokusai/boundaryParams.hpp"

namespace hokusai
{

BoundaryParams::BoundaryParams()
{
    m_radius = 0;
    m_adhesion = 0;
    m_friction = 0;
}

BoundaryParams::BoundaryParams(HReal radius, HReal adhesion, HReal friction)
{
    m_radius = radius;
    m_adhesion = adhesion;
    m_friction = friction;
    HReal eta = 0.01;
    HReal H = 0.1;
    HReal vf = std::sqrt( 2*9.81*H );
    m_soundSpeed = vf/(std::sqrt(eta));
    m_bKernel = BoundaryKernel(m_radius, m_soundSpeed);
}

BoundaryParams::~BoundaryParams()
{

}

HReal& BoundaryParams::radius()
{
    return m_radius;
}

const HReal& BoundaryParams::radius() const
{
    return m_radius;
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
    return m_bKernel;
}

const BoundaryKernel& BoundaryParams::boundaryKernel() const
{
    return m_bKernel;
}

}//namespace hokusai
