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

}//namespace hokusai
