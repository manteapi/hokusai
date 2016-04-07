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

#include "./../include/hokusai/rasterizer.hpp"

namespace hokusai
{
std::vector<Vec3i> rasterizeLine(const Vec3i& p1, const Vec3i& p2)
{
    int distance = Vec3i::max(Vec3i::abs(p1-p2)+Vec3i(1));
    std::vector<Vec3i> result(distance);

    Vec3i dP = p2-p1;
    Vec3i aP = Vec3i::abs(dP)*2;
    Vec3i sP = Vec3i::sign(dP);
    Vec3i t = p1;
    int id = 0;

    if( aP[0] >= std::max(aP[1],aP[2]) ) //x dominant
    {
        int yd = aP[1] - aP[0]/2;
        int zd = aP[2] - aP[0]/2;
        while(true)
        {
            result[id] = t;
            id++;
            if(t[0] == p2[0]) break; //end
            if(yd >= 0) //move along y
            {
                t[1] += sP[1];
                yd -= aP[0];
            }
            if(zd >= 0) // move along z
            {
                t[2] += sP[2];
                zd -= aP[0];
            }
            t[0] += sP[0]; //move along x
            yd += aP[1];
            zd += aP[2];
        }
    }
    else if( aP[1]>=std::max(aP[0],aP[2]) ) //y dominant
    {
        int xd = aP[0] - aP[1]/2;
        int zd = aP[2] - aP[1]/2;
        while(true)
        {
            result[id] = t;
            id++;
            if(t[1]==p2[1]) break; //end
            if(xd >= 0) //move along x
            {
                t[0] += sP[0];
                xd -= aP[1];
            }
            if(zd >= 0) //move along z
            {
                t[2] += sP[2];
                zd -= aP[1];
            }
            t[1] += sP[1]; //move along y
            xd += aP[0];
            zd += aP[2];
        }
    }
    else if( aP[2]>=std::max(aP[0], aP[1])) //z dominant
    {
        int xd = aP[0] - aP[2]/2;
        int yd = aP[1] - aP[2]/2;
        while(true)
        {
            result[id] = t;
            id++;
            if(t[2]==p2[2]) break;
            if(xd >= 0) // move along x
            {
                t[0] += sP[0];
                xd -= aP[2];
            }
            if(yd >= 0) // move along y
            {
                t[1] += sP[1];
                yd -= aP[2];
            }
            t[2] += sP[2]; //move along z
            xd += aP[0];
            yd += aP[1];
        }
    }

    return result;
}

void rasterizeSpan(
        const std::vector<Vec3i>& minLine,
        const std::vector<Vec3i>& maxLine,
        const int& axis,
        std::vector<Vec3i>& result, std::array< Vec3i, 2>& leftover)
{
    if(minLine[0] != maxLine[0])
    {
        std::cout << "Error" << std::endl;
        return;
    }

    result.clear();

    int minCounter=0;
    int maxCounter=0;
    while(
          (minCounter<(int)(minLine.size()))
          && (maxCounter<(int)(maxLine.size()))
          )
    {
        Vec3i p1 = minLine[minCounter];
        Vec3i p2 = maxLine[maxCounter];
        std::vector<Vec3i> t = rasterizeLine(p1, p2);

        do
        {
            maxCounter++;
        }
        while(
              (maxLine[maxCounter-1][axis] == maxLine[maxCounter][axis])
              && ((maxCounter)<(int)(maxLine.size()))
              );

        do
        {
            minCounter++;
        }
        while(
              (minLine[minCounter-1][axis] == minLine[minCounter][axis])
              && ((minCounter)<(int)(minLine.size()))
              );

        for(size_t k=0; k<t.size(); ++k)
            result.push_back(t[k]);
    }

    leftover[1] = minLine[minCounter-1];
    leftover[0] = maxLine[maxCounter-1];
}

std::vector<Vec3i> rasterizeTriangle(const Vec3i& a, const Vec3i& b, const Vec3i &c)
{
    std::vector<Vec3i> ab, ac, span, result;
    ab = rasterizeLine(a,b);
    ac = rasterizeLine(a,c);

    Vec3i dAB, dAC;
    std::array<Vec3i,2> left;
    HReal maxDistAB, maxDistAC;
    dAB = Vec3i::abs(b-a);
    maxDistAB = Vec3i::max(dAB);
    dAC = Vec3i::abs(c-a);
    maxDistAC = Vec3i::max(dAC);

    int axis = 0;

    int maxSize = std::max( maxDistAB, maxDistAC );
    if(maxSize == maxDistAB)
    {
        HReal maxAxis = maxDistAB;
        if(dAB[0] == maxAxis)
            axis = 0;
        else if(dAB[1] == maxAxis)
            axis = 1;
        else if(dAB[2] == maxAxis)
            axis = 2;
    }
    else
    {
        HReal maxAxis = maxDistAC;
        if(dAC[0] == maxAxis)
            axis = 0;
        else if(dAC[1] == maxAxis)
            axis = 1;
        else if(dAC[2] == maxAxis)
            axis = 2;
    }

    rasterizeSpan(ac, ab, axis, span, left);
    for(size_t i=0; i<span.size(); ++i)
    {
        result.push_back(span[i]);
    }
    std::vector<Vec3i> r0,r1;

    if(dAB[axis]<dAC[axis])
    {
        r0 = rasterizeLine(c,left[0]);
        r1 = rasterizeLine(c,left[1]);
    }
    else
    {
        r0 = rasterizeLine(b,left[0]);
        r1 = rasterizeLine(b,left[1]);
    }
    rasterizeSpan(r0, r1, axis, span, left);
    for(size_t i=0; i<span.size(); ++i)
    {
        result.push_back(span[i]);
    }

    return result;
}

std::vector<Vec3i> rasterizeTriangleWithoutSpan(const Vec3i& a, const Vec3i& b, const Vec3i &c)
{
    std::vector<Vec3i> ab, ac, tmp, result;
    ab = rasterizeLine(a,b);
    for(size_t i=0; i<ab.size(); ++i)
    {
        tmp = rasterizeLine(c,ab[i]);
        for(size_t j=0; j<tmp.size(); ++j)
        {
            result.push_back(tmp[j]);
        }
    }
    return result;
}

std::vector<Vec3i> overRasterizeTriangle(const Vec3i& a, const Vec3i& b, const Vec3i &c)
{
    std::vector<Vec3i> result, tmp;

    tmp = rasterizeTriangle(a,b,c);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(tmp[i]);

    tmp = rasterizeTriangle(c,a,b);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(tmp[i]);

    tmp = rasterizeTriangle(b,c,a);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(tmp[i]);

    return result;
}

std::vector<Vec3i> overRasterizeTriangleWithoutSpan(const Vec3i& a, const Vec3i& b, const Vec3i &c)
{
    std::vector<Vec3i> result, tmp;
    tmp = rasterizeTriangleWithoutSpan(a,b,c);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(tmp[i]);

    tmp = rasterizeTriangleWithoutSpan(c,a,b);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(tmp[i]);

    tmp = rasterizeTriangleWithoutSpan(b,c,a);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(tmp[i]);

    return result;
}

std::vector<Vec3i> rasterizeTriangles( const std::vector< std::array<Vec3i,3> >& triangles )
{
    std::vector<Vec3i> result;
    std::vector<Vec3i> tmpResult;
    for(size_t i=0; i<triangles.size(); ++i)
    {
        tmpResult = rasterizeTriangle( triangles[i][0], triangles[i][1], triangles[i][2] );
        for(size_t j=0; j<tmpResult.size(); ++j)
        {
            result.push_back(tmpResult[j]);
        }
    }
    return result;
}

std::vector<Vec3i> rasterizeTrianglesWithoutSpan( const std::vector< std::array<Vec3i,3> >& triangles )
{
    std::vector<Vec3i> result;
    std::vector<Vec3i> tmpResult;
    for(size_t i=0; i<triangles.size(); ++i)
    {
        tmpResult = rasterizeTriangleWithoutSpan( triangles[i][0], triangles[i][1], triangles[i][2] );
        for(size_t j=0; j<tmpResult.size(); ++j)
        {
            result.push_back(tmpResult[j]);
        }
    }
    return result;
}

std::vector<Vec3i> overRasterizeTriangles( const std::vector< std::array<Vec3i,3> >& triangles )
{
    std::vector<Vec3i> result;
    std::vector<Vec3i> tmpResult;
    for(size_t i=0; i<triangles.size(); ++i)
    {
        tmpResult.clear();
        tmpResult = overRasterizeTriangle( triangles[i][0], triangles[i][1], triangles[i][2] );
        for(size_t j=0; j<tmpResult.size(); ++j)
        {
            result.push_back(tmpResult[j]);
        }
    }
    return result;
}

std::vector<Vec3i> overRasterizeTrianglesWithoutSpan( const std::vector< std::array<Vec3i,3> >& triangles )
{
    std::vector<Vec3i> result;
    std::vector<Vec3i> tmpResult;
    for(size_t i=0; i<triangles.size(); ++i)
    {
        tmpResult.clear();
        tmpResult = overRasterizeTriangleWithoutSpan( triangles[i][0], triangles[i][1], triangles[i][2] );
        for(size_t j=0; j<tmpResult.size(); ++j)
        {
            result.push_back(tmpResult[j]);
        }
    }
    return result;
}
}
