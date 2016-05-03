#include "./../include/hokusai/fluidParams.hpp"

namespace hokusai
{

FluidParams::FluidParams()
{
    m_r = 0.0;
    m_soundSpeed = 0.0;
    m_cohesion = 0.0;
    m_viscosity = 0.0;
    m_restDensity = 0.0;
    m_h = 0.0;
    m_mass = 0.0;
}

FluidParams::FluidParams(int particleNumber, HReal volume, HReal density, HReal viscosity, HReal fluidCohesion)
{
    m_viscosity = viscosity;
    m_cohesion = fluidCohesion;
    m_restDensity = density;
    m_mass = (m_restDensity * volume) / particleNumber;
    int m_particlePerCell = 33.8; //better
    m_h = 0.5*std::pow( HReal(3.0*volume*m_particlePerCell) / HReal(4*M_PI*particleNumber), 1.0/3.0);
    m_r = std::pow( 3.0*volume/(4*M_PI*particleNumber), 1.0/3.0 );
    HReal eta = 0.01;
    HReal H = 0.1;
    HReal vf = std::sqrt( 2*9.81*H );
    m_soundSpeed = vf/(std::sqrt(eta));

}

FluidParams::~FluidParams()
{

}

HReal& FluidParams::smoothingRadius()
{
    return m_h;
}

const HReal& FluidParams::smoothingRadius() const
{
    return m_h;
}

HReal& FluidParams::radius()
{
    return m_r;
}

const HReal& FluidParams::radius() const
{
    return m_r;
}

HReal& FluidParams::mass()
{
    return m_mass;
}

const HReal& FluidParams::mass() const
{
    return m_mass;
}

HReal& FluidParams::density()
{
    return m_restDensity;
}

const HReal& FluidParams::density() const
{
    return m_restDensity;
}

HReal& FluidParams::viscosity()
{
    return m_viscosity;
}

const HReal& FluidParams::viscosity() const
{
    return m_viscosity;
}

HReal& FluidParams::cohesion()
{
    return m_cohesion;
}

const HReal& FluidParams::cohesion() const
{
    return m_cohesion;
}

HReal& FluidParams::soundSpeed()
{
    return m_soundSpeed;
}

const HReal& FluidParams::soundSpeed() const
{
    return m_soundSpeed;
}

}//namespace hokusai
