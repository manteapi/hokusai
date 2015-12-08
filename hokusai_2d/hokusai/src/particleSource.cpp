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

ParticleSource::ParticleSource(const double& startTime_, const double& endTime_, const double& delay_, const double &spacing_,
                               const Vec2d& position_, const double &orientation_, const double &scale_, const double &velocity_)
{
    position = position_;
    orientation = orientation_;
    scale = scale_;

    velocity = velocity_;

    startTime = startTime_;
    endTime = endTime_;
    delay = delay_;
    spacing = spacing_;

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
    spacing = source.spacing;

    lastTime = source.lastTime;

    init();
}

ParticleSource::ParticleSource()
{
    position = Vec2d(0.0,0.0);
    orientation = 0.0;
    scale = 1.0;

    velocity = 0.0;
    startTime = 0;
    endTime = 0;
    delay = 0;
    spacing = 0;

    lastTime = 0.0;

    init();
}

void ParticleSource::init()
{
    //Create a transformation based on the given position, radius and orientation.
//    typedef Eigen::Affine2d Affine;
//    Affine transformation;
//    transformation = Eigen::Translation2d(position[0], position[1])*Eigen::Rotation2D(angle);

    for(double dx=-scale; dx<=scale; dx+=spacing)
    {
//            Vec3d tmp3d_x = transformation*Vec3d(x, 0.0, 0.0);
            Vec2d px = position+Eigen::Rotation2Dd(orientation)*Vec2d(dx,0.0);
//            Vec3d tmp3d_v(velocity[0], velocity[1], 0.0);
//            tmp3d_v = Eigen::AngleAxisd(orientation[0], Vec3d::UnitX()) * Eigen::AngleAxisd(orientation[1],  Vec3d::UnitY())*tmp3d_v;
//            Vec2d tmp_v(tmp3d_v[0], tmp3d_v[1]);
            Particle p;
            p.x = px;
            p.v = Eigen::Rotation2Dd(orientation)*Vec2d(0,-velocity);
            p_stencil.push_back(p);
    }
}


ParticleSource::~ParticleSource()
{}

std::vector<Particle> ParticleSource::apply(const double time)
{
    int isOk = (int)((time-lastTime)/delay);
    if(isOk>0 && time<=endTime && time>=startTime)
    {
        lastTime = time;
        return p_stencil;
    }
    else
        return std::vector<Particle>();
}

}
#endif
