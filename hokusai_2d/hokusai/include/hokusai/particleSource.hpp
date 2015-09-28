#ifndef HOKUSAI_PARTICLE_SOURCE_HPP
#define HOKUSAI_PARTICLE_SOURCE_HPP

#include "utility.hpp"
#include "particle.hpp"
#include <Eigen/Geometry>

namespace hokusai
{

class ParticleSource
{
public :
    ParticleSource(const double& startTime_, const double& endTime_, const double& delay_, const double& spacing_, const Vec2d& position_=Vec2d(0,0,0), const Vec2d& orientation_=Vec2d(0,0,0), const Vec2d& scale_=Vec2d(1,1,1), const Vec2d& velocity_=Vec2d(0,0,0));
    ParticleSource(const ParticleSource& source);
    ParticleSource();
    ~ParticleSource();

    void init();

    typedef Eigen::Affine3d Affine;

    Vec2d position;
    Vec2d orientation;
    Vec2d scale;

    Vec2d velocity;
    double startTime;

    double spacing;
    double endTime;
    double delay;
    double lastTime;

    std::vector<Particle> p_stencil;

    std::vector<Particle> apply(const double time);
};

}

#endif
