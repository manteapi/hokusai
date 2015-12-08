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
{
    m_minBB = Vec2d(0.0,0.0);
    m_maxBB = Vec2d(0.0,0.0);
}

ParticleSink::ParticleSink(const Vec2d& minBB, const Vec2d& maxBB)
{
    m_minBB = minBB;
    m_maxBB = maxBB;
}

ParticleSink::ParticleSink(const ParticleSink& sink)
{
    m_minBB = sink.m_minBB;
    m_maxBB = sink.m_maxBB;
}

ParticleSink::~ParticleSink()
{}

bool ParticleSink::contain(const Vec2d& x)
{
    return ( (x[0]>=m_minBB[0] && x[0]<=m_maxBB[0]) && (x[1]>=m_minBB[1] && x[1]<=m_maxBB[1]) );
}

}
#endif
