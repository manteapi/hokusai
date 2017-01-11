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
    Mat4r m, t, rx, ry, rz;
    t.translation(position);
    rx.rotation(Vec3r(1.0,0.0,0.0), orientation[0]);
    ry.rotation(Vec3r(0.0,1.0,0.0), orientation[1]);
    rz.rotation(Vec3r(0.0,0.0,1.0), orientation[2]);
    m = t * rx * ry * rz;

    //Create a stencil of particles for the unit disc in a reference frame and orientate it using the transformation
    HReal maxRadius = Vec3r::max(scale);

    for(HReal x=-maxRadius; x<=maxRadius; x+=spacing)
    {
        for(HReal y=-maxRadius; y<=maxRadius; y+=spacing)
        {
            if( (x>=-scale[0] && x<=scale[0]) &&
                    (y>=-scale[1] && x<=scale[1]) )
            {
                Vec3r tmp_x = m*Vec3r(x, y, 0.0);
                Vec3r tmp_v = rx * ry * rz * velocity;
                Particle tmp_p;
                tmp_p.x = tmp_x;
                tmp_p.v = tmp_v;
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
