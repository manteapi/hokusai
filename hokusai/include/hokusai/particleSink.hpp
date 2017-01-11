#ifndef HOKUSAI_PARTICLE_SINK_HPP
#define HOKUSAI_PARTICLE_SINK_HPP

#include "common.hpp"
#include "gridUtility.hpp"

namespace hokusai
{

class ParticleSink
{
public :
    ~ParticleSink();
    ParticleSink();
    ParticleSink(const HReal& startTime, const HReal& endTime, const HReal& radius, const Vec3r& position);
    std::vector<int> apply(const GridUtility& gridInfo, const HReal& time);
private:
    HReal m_startTime;
    HReal m_endTime;
    HReal m_radius;
    Vec3r m_position;
};

}

#endif
