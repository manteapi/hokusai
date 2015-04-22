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

#ifndef HOKUSAI_PARTICLE_SOURCE_CPP
#define HOKUSAI_PARTICLE_SOURCE_CPP

#include "./../include/hokusai/particleSource.hpp"

namespace hokusai
{

typedef Eigen::Vector3d EigenVec;

ParticleSource::ParticleSource(const SReal& startTime_, const SReal& endTime_, const SReal& delay_, const Vec3r& position_, const Vec3r& orientation_, const Vec3r& scale_, const Vec3r& velocity_)
{
    position = position_;
    orientation = orientation_;
    scale = scale_;

    velocity = velocity_;

    startTime = startTime_;
    endTime = endTime_;
    delay = delay_;

    lastTime = 0.0;

    init();
}

ParticleSource::ParticleSource(const ParticleSource& source)
{
    position = source.position;
    orientation = source.orientation;
    scale = source.scale;

    velocity = source.velocity;
    startTime = source.startTime;
    endTime = source.endTime;
    delay = source.delay;
    lastTime = source.lastTime;

    init();
}

ParticleSource::ParticleSource()
{
    position = Vec3r(0,0,0);
    orientation = Vec3r(0,0,0);
    scale = Vec3r(1,1,1);

    velocity = Vec3r(0,0,0);
    startTime = 0;
    endTime = 0;
    delay = 0;

    lastTime = 0.0;

    init();
}

void ParticleSource::init()
{
    //Create a transformation based on the given position, radius and orientation.
    Affine transformation;
    transformation = Eigen::AngleAxisd(orientation[0], EigenVec::UnitX()) * Eigen::AngleAxisd(orientation[1],  EigenVec::UnitY()) * Eigen::AngleAxisd(orientation[2], EigenVec::UnitZ());
    transformation *= Eigen::Translation3d(position[0], position[1], position[2]);
    transformation *= Eigen::Scaling(scale[0], scale[1], scale[2]);

    //Create a stencil of particles for the unit disc in a reference frame and orientate it using the transformation
    p_stencil.resize(4);
    EigenVec tmp_x;
    tmp_x = transformation*EigenVec(0,0,0);
    p_stencil[0].x = Vec3r(tmp_x[0], tmp_x[1], tmp_x[2]);
    p_stencil[0].v = velocity;

    tmp_x = transformation*EigenVec(1,0,0);
    p_stencil[1].x = Vec3r(tmp_x[0], tmp_x[1], tmp_x[2]);
    p_stencil[1].v = velocity;

    tmp_x = transformation*EigenVec(0,1,0);
    p_stencil[2].x = Vec3r(tmp_x[0], tmp_x[1], tmp_x[2]);
    p_stencil[2].v = velocity;

    tmp_x = transformation*EigenVec(1,1,0);
    p_stencil[3].x = Vec3r(tmp_x[0], tmp_x[1], tmp_x[2]);
    p_stencil[3].v = velocity;
}


ParticleSource::~ParticleSource()
{}

std::vector<Particle> ParticleSource::apply(const SReal time)
{
    int isOk = (int)((time-lastTime)/delay);
    if(isOk>0)
    {
        lastTime = time;
        return p_stencil;
    }
    else
        return std::vector<Particle>();
}

}
#endif
