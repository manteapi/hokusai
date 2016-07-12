#ifndef HOKUSAI_FLUID_PARAMS
#define HOKUSAI_FLUID_PARAMS

#include "common.hpp"

namespace hokusai
{

class FluidParams
{
public:
    FluidParams();
    FluidParams(const HReal& particleNumber, const HReal& volume, const HReal& restDensity);
    FluidParams(const FluidParams& params);
    ~FluidParams();
    HReal& mass();
    const HReal& mass() const;
    HReal& restDensity();
    const HReal& restDensity() const;
    HReal& smoothingRadius();
    const HReal& smoothingRadius() const;

private:
    HReal m_smoothingRadius;
    HReal m_mass;
    HReal m_restDensity;
};

} //namespace hokusai

#endif //HOKUSAI_FLUID_PARAMS
