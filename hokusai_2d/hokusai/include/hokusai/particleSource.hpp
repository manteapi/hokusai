#ifndef HOKUSAI_PARTICLE_SOURCE_HPP
#define HOKUSAI_PARTICLE_SOURCE_HPP

#include <Eigen/Geometry>
#include "particle.hpp"
#include "utility.hpp"

namespace hokusai
{

class ParticleSource
{
public :
    ParticleSource(const double& startTime_, const double& endTime_, const double& delay_, const double& spacing_,
                   const Vec2d& position_=Vec2d(0.0,0.0), const double& orientation_=0.0,
                   const double& scale_= 1.0, const double& velocity_= 0.0);
    ParticleSource(const ParticleSource& source);
    ParticleSource();
    ~ParticleSource();

    void init();

    Vec2d position;
    double orientation;
    double scale;

    double velocity;
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
