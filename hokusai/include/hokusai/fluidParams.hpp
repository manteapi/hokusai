#ifndef HOKUSAI_FLUID_PARAMS
#define HOKUSAI_FLUID_PARAMS

#include "common.hpp"
#include "kernel.hpp"

namespace hokusai
{

class FluidParams
{
public:
    FluidParams();
    FluidParams(const HReal& particleNumber, const HReal& volume, const HReal& restDensity, const HReal &viscosity,
                const HReal& cohesion);
    FluidParams(const FluidParams& params);
    ~FluidParams();
    HReal& mass();
    const HReal& mass() const;
    HReal& restDensity();
    const HReal& restDensity() const;
    HReal& smoothingRadius();
    const HReal& smoothingRadius() const;
    HReal& viscosity();
    const HReal& viscosity() const;
    HReal& soundSpeed();
    const HReal& soundSpeed() const;
    HReal& cohesion();
    const HReal& cohesion() const;
    MonaghanKernel& monaghanKernel();
    const MonaghanKernel& monaghanKernel() const;
    AkinciKernel& akinciKernel();
    const AkinciKernel& akinciKernel() const;
private:
    HReal m_viscosity;
    HReal m_smoothingRadius;
    HReal m_mass;
    HReal m_restDensity;
    HReal m_soundSpeed;
    HReal m_cohesion;
    MonaghanKernel m_pKernel;
    AkinciKernel m_aKernel;
};

} //namespace hokusai

#endif //HOKUSAI_FLUID_PARAMS
