#include "./../include/hokusai/fluidParams.hpp"

namespace hokusai
{

FluidParams::FluidParams()
{}

FluidParams::FluidParams(const HReal& particleNumber, const HReal& volume, const HReal& restDensity,
                         const HReal& viscosity, const HReal &cohesion)
{
    m_cohesion = cohesion;
    m_viscosity = viscosity;
    m_restDensity = restDensity;
    m_mass = (m_restDensity * volume) / particleNumber;
    HReal particlePerCell = 33.8; //better
    m_smoothingRadius = 0.5*pow( HReal(3*volume*particlePerCell) / HReal(4*M_PI*particleNumber), 1.0/3.0);

    HReal eta = 0.01;
    HReal H = 0.1;
    HReal vf = sqrt( 2*9.81*H );
    m_soundSpeed = vf/(sqrt(eta));

    m_pKernel = MonaghanKernel(m_smoothingRadius);
    m_aKernel = AkinciKernel(2.0*m_smoothingRadius);
}

FluidParams::FluidParams(const FluidParams& params)
{
    m_viscosity = params.viscosity();
    m_cohesion = params.cohesion();
    m_restDensity = params.restDensity();
    m_mass = params.mass();
    m_smoothingRadius = params.mass();
    m_pKernel = params.monaghanKernel();
    m_aKernel = params.akinciKernel();
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

const HReal& FluidParams::viscosity() const
{
    return m_viscosity;
}
HReal& FluidParams::viscosity()
{
    return m_viscosity;
}

HReal& FluidParams::soundSpeed()
{
    return m_soundSpeed;
}

const HReal& FluidParams::soundSpeed() const
{
    return m_soundSpeed;
}

HReal& FluidParams::cohesion()
{
    return m_cohesion;
}

const HReal& FluidParams::cohesion() const
{
    return m_cohesion;
}

MonaghanKernel& FluidParams::monaghanKernel()
{
    return m_pKernel;
}

const MonaghanKernel& FluidParams::monaghanKernel() const
{
    return m_pKernel;
}

AkinciKernel& FluidParams::akinciKernel()
{
    return m_aKernel;
}

const AkinciKernel& FluidParams::akinciKernel() const
{
    return m_aKernel;
}

} // namespace hokusai
