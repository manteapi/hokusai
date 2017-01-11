#include "./../include/hokusai/solverParams.hpp"

namespace hokusai
{

SolverParams::~SolverParams(){}

SolverParams::SolverParams(){}

SolverParams::SolverParams(const HReal& timeStep, const int &maxPressureSolveIterationNb,
             const HReal& maxDensityError)
{
    m_timeStep = timeStep;
    m_maxPressureSolveIterationNb = maxPressureSolveIterationNb;
    m_maxDensityError = maxDensityError;
    m_averageDensity = 0.0;
}

SolverParams::SolverParams(const SolverParams& params)
{
    m_timeStep = params.timeStep();
    m_averageDensity = params.averageDensity();
    m_maxPressureSolveIterationNb = params.maxPressureSolveIterationNb();
    m_maxDensityError = params.maxDensityError();
}

HReal& SolverParams::timeStep()
{
    return m_timeStep;
}

const HReal& SolverParams::timeStep() const
{
    return m_timeStep;
}

HReal& SolverParams::averageDensity()
{
    return m_averageDensity;
}

const HReal& SolverParams::averageDensity() const
{
    return m_averageDensity;
}

int &SolverParams::maxPressureSolveIterationNb()
{
    return m_maxPressureSolveIterationNb;
}

const int& SolverParams::maxPressureSolveIterationNb() const
{
    return m_maxPressureSolveIterationNb;
}

HReal& SolverParams::maxDensityError()
{
    return m_maxDensityError;
}

const HReal& SolverParams::maxDensityError() const
{
    return m_maxDensityError;
}

} //namespace hokusai
