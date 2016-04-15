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

#ifndef HOKUSAI_UTILS_INL
#define HOKUSAI_UTILS_INL

#include "utils.hpp"
#include <vector>
#include "write_bmp.hpp"
#include "particle.hpp"
#include "boundary.hpp"
#include "common.hpp"
#include <fstream>
#include <sstream>

namespace hokusai
{

template<typename ParticleT>
void write_frame(ParticleContainer<ParticleT>& particles, int step, HReal offset)
{
    int width = 1024;
    int height = 700;
    HReal winrate = height/(HReal)width;

    static HReal *image = new HReal[width*height*4];
    static unsigned char *buffer = new unsigned char[width*height*4];
    for( int i=0; i<width; i++ ) {
        for( int j=0; j<height; j++ ) {
            for( int c=0; c<3; c++ ) {
                image[4*(i+j*width)+c] = 0.0; //black
            }
            image[4*(i+j*width)+3] = 1.0; //opacity
        }
    }

    HReal DENSITY = 0.5;
    int N = 32;
    HReal blueColor[] = { 0.3, 0.5, 0.8, exp(-DENSITY*N/23.0f) };

    static HReal R[3][3];
    bool firstTime = true;
    if( firstTime ) {
        //buildRotationMatrix( -0.2, 0.2, R );
        buildRotationMatrix( 0, 0, R );
        firstTime = false;
    }

    // Simple Point-based Rendering
    //HReal eye = 2.0;
    //HReal offset = 0.3;
    HReal eye = 2.0;

    for( size_t n=0; n<particles.size(); n++ )
    {
        //if( particles[n]->type == FLUID ) %{
        //HReal p[3] = { particles[n]->p[0]-0.5, particles[n]->p[1]-0.5, particles[n]->p[2]-0.5 };
        HReal p[3] = { (HReal)(particles[n].x[0]-0.5), (HReal)(particles[n].x[1]-0.5), (HReal)(particles[n].x[2]-0.5) };
        transform(p,R);
        HReal z = offset + 0.5 + p[2];
        HReal x = eye/(eye+z)*(p[0]-0.4);
        HReal y = eye/(eye+z)*(p[1]+0.25);
        if( x > -0.5 && x < 0.5 && y > -winrate/2.0 && y < winrate/2.0 ) {
            int i = width*(x+0.5);
            int j = width*(y+winrate/2.0);
            int w = 1;
            for( int ni=i-w; ni<i+1; ni++ ) {
                for( int nj=j-w; nj<j+1; nj++ ) {
                    if( ni>0 && ni<width && nj>0 && nj<height ) {
                        for( int c=0; c<3; c++ ) {
                            image[4*(ni+nj*width)+c] = blueColor[3]*blueColor[c] + (1.0-blueColor[3])*image[4*(ni+nj*width)+c];
                        }
                    }
                }
            }
        }
        //}
    }

    for( int n=0; n<4*width*height; n++ ) {
        buffer[n] = 255.0*fmin(1.0,image[n]);
    }

    char name[256];
    sprintf( name, "frame_%d.bmp", step );
    write_bmp( name, buffer, width, height, true );
}

}

#endif //HOKUSAI_UTILS_INL
