#ifndef HOKUSAI_PARTICLE_SOURCE_HPP
#define HOKUSAI_PARTICLE_SOURCE_HPP

#include "Vec.hpp"

namespace hokusai
{

class ParticleSource
{
public :
    ParticleSource();
    ~ParticleSource();

    Vec3<double> position;
    Vec3<double> velocity;
    double radius, startTime, endTime;
};

}

#endif
