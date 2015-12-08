#ifndef HOKUSAI_PARTICLE_SINK_HPP
#define HOKUSAI_PARTICLE_SINK_HPP

#include "utility.hpp"

namespace hokusai
{

class ParticleSink
{
public :
    ParticleSink();
    ParticleSink(const Vec2d& minBB, const Vec2d& maxBB);
    ParticleSink(const ParticleSink& sink);
    ~ParticleSink();
    Vec2d m_minBB, m_maxBB;
    bool contain(const Vec2d& x);
};

}

#endif
