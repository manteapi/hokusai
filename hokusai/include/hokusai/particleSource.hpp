#ifndef HOKUSAI_PARTICLE_SOURCE_HPP
#define HOKUSAI_PARTICLE_SOURCE_HPP

#include <aljabr/Vec.hpp>

namespace hokusai
{

    typedef aljabr::Vec3<double> Vec3r;

class ParticleSource
{
public :
    ParticleSource();
    ~ParticleSource();

    Vec3r position;
    Vec3r velocity;
    double radius, startTime, endTime;
};

}

#endif
