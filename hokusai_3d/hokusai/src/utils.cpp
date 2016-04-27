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

#include "../include/hokusai/utils.hpp"
#include <bitset>

namespace hokusai
{

int mortonNumber( std::array<int,3>& index )
{
    std::bitset<64> bitsMorton;
    std::array< std::bitset<64>, 3> bitsIndex;
    for(size_t i=0; i<bitsIndex.size(); ++i) bitsIndex[i] = std::bitset<64>(index[i]);
    for(size_t i=0; i<bitsMorton.size()/3; ++i)
        for(size_t j=0; j<bitsIndex.size(); ++j)
            bitsMorton[3*i+j] = bitsIndex[j][i];
    int result = bitsMorton.to_ullong();
    return result;
}

void write(const char * filename, std::vector<HReal> data)
{
    std::ofstream outputFile;
    outputFile.open(filename);
    outputFile.precision(16);
    for(unsigned int i=0; i <data.size(); ++i)
    {
        outputFile << data[i] <<"\n";
    }
    outputFile.close();
}

void write(const char * filename, std::vector<Vec3r > data)
{
    std::ofstream outputFile;
    outputFile.open(filename);
    outputFile.precision(16);
    for(unsigned int i=0; i <data.size(); ++i)
    {
        outputFile << data[i][0] << " " << data[i][1] <<" " <<  data[i][2] << "\n";
    }
    outputFile.close();
}

static void buildRotationMatrix( HReal xrad, HReal yrad, HReal R[3][3] ) {
    HReal Rx[3][3] = { {1,0,0},{0,cos(xrad),-sin(xrad)},{0,sin(xrad),cos(xrad)} };
    HReal Ry[3][3] = { {cos(yrad),0,sin(yrad)}, {0,1,0}, {-sin(yrad),0,cos(yrad)} };
    HReal Rtmp[3][3];
    for( int i=0; i<3; i++ ) {
        for( int j=0; j<3; j++ ) {
            R[i][j] = Rtmp[i][j] = Rx[i][j];
        }
    }

    for( int i=0; i<3; i++ ) {
        for( int j=0; j<3; j++ ) {
            R[i][j] = 0.0;
            for( int k=0; k<3; k++ ) R[i][j] += Rtmp[i][k]*Ry[k][j];
        }
    }
}

static void transform( HReal p[3], HReal R[3][3] ) {
    HReal p0[3] = { p[0], p[1], p[2] };
    for( int i=0; i<3; i++ ) {
        p[i] = 0.0;
        for( int k=0; k<3; k++ ) {
            p[i] += R[i][k]*p0[k];
        }
    }
}

void write_frame(std::vector<Particle>& particles, int step, HReal offset)
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

}//namespace hokusai
