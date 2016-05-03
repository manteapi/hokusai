#ifndef HOKUSAI_FLUID_PARAMS_HPP
#define HOKUSAI_FLUID_PARAMS_HPP

#include "common.hpp"

namespace hokusai
{

class FluidParams
{
public:
    FluidParams();
    FluidParams(int particleNumber, HReal volume, HReal density,
                HReal viscosity=0.1, HReal fluidCohesion=0.05);
    ~FluidParams();

    HReal& smoothingRadius();
    const HReal& smoothingRadius() const;

    HReal& radius();
    const HReal& radius() const;

    HReal& mass();
    const HReal& mass() const;

    HReal& density();
    const HReal& density() const;

    HReal& viscosity();
    const HReal& viscosity() const;

    HReal& cohesion();
    const HReal& cohesion() const;

    HReal& soundSpeed();
    const HReal& soundSpeed() const;

private:
    HReal m_soundSpeed;
    HReal m_cohesion;
    HReal m_viscosity;
    HReal m_restDensity;
    HReal m_h;
    HReal m_r;
    HReal m_mass;
};

}//namespace hokusai

#endif//HOKUSAI_FLUID_PARAMS_HPP
