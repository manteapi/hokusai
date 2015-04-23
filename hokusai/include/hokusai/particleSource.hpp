#ifndef HOKUSAI_PARTICLE_SOURCE_HPP
#define HOKUSAI_PARTICLE_SOURCE_HPP

#include <aljabr/Vec.hpp>
#include <Eigen/Geometry>
#include "particle.hpp"

namespace hokusai
{

typedef aljabr::Vec3<double> Vec3r;
typedef double SReal;
class ParticleSource
{
public :
    ParticleSource(const SReal& startTime_, const SReal& endTime_, const SReal& delay_, const SReal& spacing_, const Vec3r& position_=Vec3r(0,0,0), const Vec3r& orientation_=Vec3r(0,0,0), const Vec3r& scale_=Vec3r(1,1,1), const Vec3r& velocity_=Vec3r(0,0,0));
    ParticleSource(const ParticleSource& source);
    ParticleSource();
    ~ParticleSource();

    void init();

    typedef Eigen::Affine3d Affine;

    Vec3r position;
    Vec3r orientation;
    Vec3r scale;

    Vec3r velocity;
    SReal startTime;

    SReal spacing;
    SReal endTime;
    SReal delay;
    SReal lastTime;

    std::vector<Particle> p_stencil;

    std::vector<Particle> apply(const SReal time);
};

}

#endif
