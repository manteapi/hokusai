#ifndef HOKUSAI_IISPH_HPP
#define HOKUSAI_IISPH_HPP

#include "kernel.hpp"

namespace hokusai
{

class IISPHParams
{
public:
    IISPHParams(const HReal& fluidSmoothingRadius, const HReal& boundarySmoothingRadius, const HReal &soundSpeed,
                const HReal& maxDensityError, const int& maxPressureSolverIterationNb);
    IISPHParams();
    IISPHParams(const IISPHParams& params);
    ~IISPHParams();

    HReal& averageDensity();
    const HReal& averageDensity() const;

    HReal& maxDensityError();
    const HReal& maxDensityError() const;

    int& maxPressureSolveIterationNb();
    const int& maxPressureSolveIterationNb() const;

    AkinciKernel& akinciKernel();
    const AkinciKernel& akinciKernel() const;

    MonaghanKernel& monaghanKernel();
    const MonaghanKernel& monaghanKernel() const;

    BoundaryKernel& boundaryKernel();
    const BoundaryKernel& boundaryKernel() const;

private:
    HReal m_rhoAvg;                     /*!< Average density at the last pressure solve iteration */
    HReal m_maxRhoError;                /*!< Maximum density error */
    int m_maxPressureIterNb;
    AkinciKernel m_aKernel;
    MonaghanKernel m_pKernel;
    BoundaryKernel m_bKernel;
};

}//namespace hokusai

#endif //HOKUSAI_IISPH_HPP
