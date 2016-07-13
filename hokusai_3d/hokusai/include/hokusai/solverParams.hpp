#ifndef HOKUSAI_SOLVER_PARAMS
#define HOKUSAI_SOLVER_PARAMS

#include "common.hpp"

namespace hokusai
{

class SolverParams
{
public:
    SolverParams();
    SolverParams(const HReal& timeStep, const int& maxPressureSolveIterationNb,
                 const HReal& maxDensityError);
    SolverParams(const SolverParams& params);
    ~SolverParams();
    HReal& timeStep();
    const HReal& timeStep() const;
    HReal& averageDensity();
    const HReal& averageDensity() const;
    int& maxPressureSolveIterationNb();
    const int &maxPressureSolveIterationNb() const;
    HReal& maxDensityError();
    const HReal& maxDensityError() const;
private:
    HReal m_timeStep;
    HReal m_averageDensity;
    int m_maxPressureSolveIterationNb;
    HReal m_maxDensityError;
};

} //namespace hokusai

#endif //HOKUSAI_SOLVER_PARAMS
