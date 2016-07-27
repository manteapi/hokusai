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

void buildRotationMatrix( HReal xrad, HReal yrad, HReal R[3][3] )
{
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

void transform( HReal p[3], HReal R[3][3] )
{
    HReal p0[3] = { p[0], p[1], p[2] };
    for( int i=0; i<3; i++ ) {
        p[i] = 0.0;
        for( int k=0; k<3; k++ ) {
            p[i] += R[i][k]*p0[k];
        }
    }
}

}//namespace hokusai
