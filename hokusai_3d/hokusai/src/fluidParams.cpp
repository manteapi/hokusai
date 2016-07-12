#include "./../include/hokusai/fluidParams.hpp"

namespace hokusai
{

FluidParams::FluidParams()
{

}

FluidParams::FluidParams(const HReal& particleNumber, const HReal& volume, const HReal& restDensity)
{
    m_restDensity = restDensity;
    m_mass = (m_restDensity * volume) / particleNumber;
    HReal particlePerCell = 33.8; //better
    m_smoothingRadius = 0.5*pow( HReal(3*volume*particlePerCell) / HReal(4*M_PI*particleNumber), 1.0/3.0);
}

FluidParams::FluidParams(const FluidParams& params)
{
    m_restDensity = params.restDensity();
    m_mass = params.mass();
    m_smoothingRadius = params.mass();
}

FluidParams::~FluidParams()
{

}

HReal& FluidParams::mass()
{
    return m_mass;
}

const HReal& FluidParams::mass() const
{
    return m_mass;
}

HReal& FluidParams::restDensity()
{
    return m_restDensity;
}

const HReal& FluidParams::restDensity() const
{
    return m_restDensity;
}

HReal& FluidParams::smoothingRadius()
{
    return m_smoothingRadius;
}

const HReal& FluidParams::smoothingRadius() const
{
    return m_smoothingRadius;
}

} // namespace hokusai
