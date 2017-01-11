/*
* Copyright 2015 Pierre-Luc Manteaux
*
*This file is part of Hokusai.
*
*Hokusai is free software: you can redistribute it and/or modify
*it under the terms of the GNU General Public License as published by
*the Free Software Foundation, either version 3 of the License, or
*(at your option) any later version.
*
*Hokusai is distributed in the hope that it will be useful,
*but WITHOUT ANY WARRANTY; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*GNU General Public License for more details.
*
*You should have received a copy of the GNU General Public License
*along with Hokusai.  If not, see <http://www.gnu.org/licenses/>.
*
* Author : Pierre-Luc Manteaux
* Contact : pierre-luc.manteaux@inria.fr
*/

#ifndef HOKUSAI_PARTICLE_SINK_CPP
#define HOKUSAI_PARTICLE_SINK_CPP

#include "./../include/hokusai/particleSink.hpp"

namespace hokusai
{

ParticleSink::ParticleSink()
{}

ParticleSink::~ParticleSink()
{}

ParticleSink::ParticleSink(const HReal& startTime, const HReal& endTime, const HReal& radius, const Vec3r& position)
{
    m_startTime = startTime;
    m_endTime = endTime;
    m_position = position;
    m_radius = radius;
}

std::vector<int> ParticleSink::apply(const GridUtility& gridInfo, const HReal& time)
{
    std::vector<int> result;
    if(time>=m_startTime && time <=m_endTime)
    {
        gridInfo.get27Neighbors(result, m_position, m_radius);
    }
    return result;
}

}
#endif
