#include "./../include/hokusai/iisph.hpp"

namespace hokusai
{

IISPHParams::IISPHParams()
{
    m_pKernel = MonaghanKernel();
    m_aKernel = AkinciKernel();
    m_bKernel = BoundaryKernel();
}

IISPHParams::IISPHParams( const IISPHParams& params)
{
    m_pKernel = MonaghanKernel(params.monaghanKernel());
    m_aKernel = AkinciKernel(params.akinciKernel());
    m_bKernel = BoundaryKernel(params.boundaryKernel());
}

IISPHParams::IISPHParams(const HReal& fluidSmoothingRadius, const HReal& boundarySmoothingRadius, const HReal& soundSpeed, const HReal &maxDensityError, const int &maxPressureSolverIterationNb)
{
    m_maxPressureIterNb = maxPressureSolverIterationNb;
    m_maxRhoError = maxDensityError;
    m_pKernel = MonaghanKernel(fluidSmoothingRadius);
    m_aKernel = AkinciKernel(2.0*fluidSmoothingRadius);
    m_bKernel = BoundaryKernel(boundarySmoothingRadius, soundSpeed);
}

IISPHParams::~IISPHParams()
{

}

AkinciKernel& IISPHParams::akinciKernel()
{
    return m_aKernel;
}

const AkinciKernel& IISPHParams::akinciKernel() const
{
    return m_aKernel;
}

MonaghanKernel& IISPHParams::monaghanKernel()
{
    return m_pKernel;
}

const MonaghanKernel& IISPHParams::monaghanKernel() const
{
    return m_pKernel;
}

BoundaryKernel& IISPHParams::boundaryKernel()
{
    return m_bKernel;
}

const BoundaryKernel& IISPHParams::boundaryKernel() const
{
    return m_bKernel;
}

HReal& IISPHParams::averageDensity()
{
    return m_rhoAvg;
}

const HReal& IISPHParams::averageDensity() const
{
    return m_rhoAvg;
}

HReal& IISPHParams::maxDensityError()
{
    return m_maxRhoError;
}

const HReal& IISPHParams::maxDensityError() const
{
    return m_maxRhoError;
}

int& IISPHParams::maxPressureSolveIterationNb()
{
    return m_maxPressureIterNb;
}

const int& IISPHParams::maxPressureSolveIterationNb() const
{
    return m_maxPressureIterNb;
}

}//namespace hokusai
