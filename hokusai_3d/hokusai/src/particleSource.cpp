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

ParticleSource::ParticleSource(const HReal& startTime_, const HReal& endTime_, const HReal& delay_, const HReal &spacing_, const Vec3r& position_, const Vec3r& orientation_, const Vec3r& scale_, const Vec3r& velocity_)
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
    position = Vec3r(0,0,0);
    orientation = Vec3r(0,0,0);
    scale = Vec3r(1,1,1);

    velocity = Vec3r(0,0,0);
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
    Affine transformation;
    transformation = Eigen::Translation3d(position[0], position[1], position[2]) * Eigen::AngleAxisd(orientation[0], EigenVec::UnitX()) * Eigen::AngleAxisd(orientation[1],  EigenVec::UnitY()) * Eigen::AngleAxisd(orientation[2], EigenVec::UnitZ());

    //Create a stencil of particles for the unit disc in a reference frame and orientate it using the transformation
    EigenVec radiusVec  = Eigen::Scaling(scale[0], scale[1], scale[2])*EigenVec(1.0,1.0,1.0);
    HReal maxRadius = radiusVec.maxCoeff();

    for(HReal x=-maxRadius; x<=maxRadius; x+=spacing)
    {
        for(HReal y=-maxRadius; y<=maxRadius; y+=spacing)
        {
            if( (x>=-radiusVec[0] && x<=radiusVec[0]) &&
                    (y>=-radiusVec[1] && x<=radiusVec[1]) )
            {
                EigenVec tmp_x = transformation*EigenVec(x, y, 0.0);
                EigenVec tmp_v = Eigen::AngleAxisd(orientation[0], EigenVec::UnitX()) * Eigen::AngleAxisd(orientation[1],  EigenVec::UnitY()) * Eigen::AngleAxisd(orientation[2], EigenVec::UnitZ())*EigenVec(velocity[0], velocity[1], velocity[2]);
                Particle tmp_p;
                tmp_p.x = Vec3r(tmp_x[0], tmp_x[1], tmp_x[2]);
                tmp_p.v = Vec3r(tmp_v[0], tmp_v[1], tmp_v[2]);
                p_stencil.push_back(tmp_p);
            }
        }
    }
}


ParticleSource::~ParticleSource()
{}

std::vector<Particle> ParticleSource::apply(const HReal time)
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
