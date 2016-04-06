#ifndef HOKUSAI_PARTICLE_SOURCE_HPP
#define HOKUSAI_PARTICLE_SOURCE_HPP

#include <aljabr/AljabrCore>
#include <Eigen/Geometry>
#include "particle.hpp"

namespace hokusai
{

class ParticleSource
{
public :
    ParticleSource(const HReal& startTime_, const HReal& endTime_, const HReal& delay_, const HReal& spacing_, const Vec3r& position_=Vec3r(0,0,0), const Vec3r& orientation_=Vec3r(0,0,0), const Vec3r& scale_=Vec3r(1,1,1), const Vec3r& velocity_=Vec3r(0,0,0));
    ParticleSource(const ParticleSource& source);
    ParticleSource();
    ~ParticleSource();

    void init();

    typedef Eigen::Affine3d Affine;

    Vec3r position;
    Vec3r orientation;
    Vec3r scale;

    Vec3r velocity;
    HReal startTime;

    HReal spacing;
    HReal endTime;
    HReal delay;
    HReal lastTime;

    std::vector<Particle> p_stencil;

    std::vector<Particle> apply(const HReal time);
};

}

#endif
