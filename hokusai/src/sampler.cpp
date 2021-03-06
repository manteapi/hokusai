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

#include "./../include/hokusai/sampler.hpp"
#include "./../include/hokusai/gridUtility.hpp"

namespace hokusai
{
/*
   Calculate the line segment PaPb that is the shortest route between
   two lines P1P2 and P3P4. Calculate also the values of mua and mub where
      Pa = P1 + mua (P2 - P1)
      Pb = P3 + mub (P4 - P3)
   Return FALSE if no solution exists.
*/
bool LineLineIntersect(
        const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const Vec3r& p4, Vec3r& pa, Vec3r& pb,
        HReal & mua, HReal & mub)
{
    Vec3r p13,p43,p21;
    HReal d1343,d4321,d1321,d4343,d2121;
    HReal numer,denom;

    p13[0] = p1[0] - p3[0];
    p13[1] = p1[1] - p3[1];
    p13[2] = p1[2] - p3[2];
    p43[0] = p4[0] - p3[0];
    p43[1] = p4[1] - p3[1];
    p43[2] = p4[2] - p3[2];
    if (std::abs(p43[0]) < std::numeric_limits<HReal>::epsilon() && std::abs(p43[1]) < std::numeric_limits<HReal>::epsilon() && std::abs(p43[2]) < std::numeric_limits<HReal>::epsilon())
        return false;
    p21[0] = p2[0] - p1[0];
    p21[1] = p2[1] - p1[1];
    p21[2] = p2[2] - p1[2];
    if (std::abs(p21[0]) < std::numeric_limits<HReal>::epsilon() && std::abs(p21[1]) < std::numeric_limits<HReal>::epsilon() && std::abs(p21[2]) < std::numeric_limits<HReal>::epsilon())
        return false;

    d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2];
    d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2];
    d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2];
    d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
    d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];

    denom = d2121 * d4343 - d4321 * d4321;
    if (std::abs(denom) < std::numeric_limits<HReal>::epsilon())
        return false;
    numer = d1343 * d4321 - d1321 * d4343;

    mua = numer / denom;
    mub = (d1343 + d4321 * (mua)) / d4343;

    pa[0] = p1[0] + mua * p21[0];
    pa[1] = p1[1] + mua * p21[1];
    pa[2] = p1[2] + mua * p21[2];
    pb[0] = p3[0] + mub * p43[0];
    pb[1] = p3[1] + mub * p43[1];
    pb[2] = p3[2] + mub * p43[2];

    return true;
}

bool LineIntersect(const Vec2r& p1, const Vec2r& p2, const Vec2r& p3, const Vec2r& p4, Vec2r& p)
{
    HReal mua,mub;
    HReal denom,numera,numerb;

    denom  = (p4[1]-p3[1]) * (p2[0]-p1[0]) - (p4[0]-p3[0]) * (p2[1]-p1[1]);
    numera = (p4[0]-p3[0]) * (p1[1]-p3[1]) - (p4[1]-p3[1]) * (p1[0]-p3[0]);
    numerb = (p2[0]-p1[0]) * (p1[1]-p3[1]) - (p2[1]-p1[1]) * (p1[0]-p3[0]);

    /* Are the line coincident? */
    if (std::abs(numera) < std::numeric_limits<HReal>::epsilon() && std::abs(numerb) < std::numeric_limits<HReal>::epsilon() && std::abs(denom) < std::numeric_limits<HReal>::epsilon()) {
        p[0] = (p1[0] + p2[0]) / 2;
        p[1] = (p1[1] + p2[1]) / 2;
        return true;
    }

    /* Are the line parallel */
    if (std::abs(denom) < std::numeric_limits<HReal>::epsilon()) {
        p[0] = 0;
        p[1] = 0;
        std::cerr << "Parallel Line" << std::endl;
        return false;
    }

    /* Is the intersection along the the segments */
    mua = numera / denom;
    mub = numerb / denom;
    /*if (mua < 0 || mua > 1 || mub < 0 || mub > 1) {
       std::cerr << "Intersection out of segment" << std::endl;
      p[0] = 0;
      p[1] = 0;
      return false;
   }*/
    p[0] = p1[0] + mua * (p2[0] - p1[0]);
    p[1] = p1[1] + mua * (p2[1] - p1[1]);
    return true;
}


bool AkinciEdgeSampling( const Vec3r& p1, const Vec3r& p2, const HReal& particleDiameter, std::vector< Vec3r >& samples)
{
    samples.clear();

    Vec3r edge = p2-p1;
    HReal edgesL = edge.length();
    int pNumber = std::floor(edgesL/particleDiameter);
    Vec3r pe = edge/(HReal)pNumber, p;

    for(int j=1; j<pNumber; ++j)
    {
        p = p1 + (HReal)j*pe;
        samples.push_back(p);
    }

    return true;
}

bool AkinciFullTriangleSampling( const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const HReal& particleDiameter, std::vector< Vec3r >& samples)
{
    std::vector< Vec3r > tmp;

    samples.push_back(p1);
    samples.push_back(p2);
    samples.push_back(p3);

    bool success = AkinciTriangleSampling(p1,p2,p3,particleDiameter, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        samples.push_back(tmp[i]);

    AkinciEdgeSampling(p1,p2,particleDiameter, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        samples.push_back(tmp[i]);

    AkinciEdgeSampling(p1,p3,particleDiameter, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        samples.push_back(tmp[i]);

    AkinciEdgeSampling(p2,p3,particleDiameter, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        samples.push_back(tmp[i]);
   
    return success;
}

bool AkinciTriangleSampling( const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const HReal& particleDiameter, std::vector< Vec3r >& samples)
{

    std::array< Vec3r, 3> v = {{p1, p2, p3}};
    std::array< Vec3r, 3 > edgesV = {{v[1]-v[0], v[2]-v[1], v[0]-v[2]}};
    std::array< Vec2i, 3 > edgesI = {{Vec2i(0,1), Vec2i(1,2), Vec2i(2,0)}};
    std::array< HReal, 3> edgesL = {{ edgesV[0].length(), edgesV[1].length(), edgesV[2].length() }};
    samples.clear();

    //Edges
    int pNumber=0;
    Vec3r pe(0,0,0);
    for(int i=0; i<3; ++i)
    {
        pNumber = std::floor(edgesL[i]/particleDiameter);
        pe = edgesV[i]/(HReal)pNumber;
        for(int j=0; j<pNumber; ++j)
        {
            Vec3r p = v[edgesI[i][0]] + (HReal)j*pe;
            //samples.push_back(p);
        }
    }

    //Triangles
    int sEdge=-1,lEdge=-1;
    HReal maxL = -std::numeric_limits<HReal>::max();
    HReal minL = std::numeric_limits<HReal>::max();
    for(int i=0; i<3; ++i)
    {
        if(edgesL[i]>maxL)
        {
            maxL = edgesL[i];
            lEdge = i;
        }
        if(edgesL[i]<minL)
        {
            minL = edgesL[i];
            sEdge = i;
        }
    }
    Vec3r cross, normal;
    cross = Vec3r::crossProduct(edgesV[lEdge], edgesV[sEdge]);
    normal = Vec3r::crossProduct(edgesV[sEdge], cross);
    normal.normalize();


    std::array<bool, 3> findVertex = {{true, true, true}};
    findVertex[edgesI[sEdge][0]] = false;
    findVertex[edgesI[sEdge][1]] = false;
    int thirdVertex = -1;
    for(size_t i=0; i<findVertex.size(); ++i)
        if(findVertex[i]==true)
            thirdVertex = i;
    Vec3r tmpVec  = v[thirdVertex] - v[edgesI[sEdge][0]];
    HReal sign = Vec3r::dotProduct(normal, tmpVec);
    if(sign<0)
        normal = -normal;

    HReal triangleHeight = std::abs(Vec3r::dotProduct(normal, edgesV[lEdge]));
    int sweepSteps = triangleHeight/particleDiameter;
    bool success = false;

    Vec3r sweepA, sweepB, i1, i2, o1, o2;
    HReal m1, m2;
    int edge1,edge2;
    edge1 = (sEdge+1)%3;
    edge2 = (sEdge+2)%3;

    for(int i=1; i<sweepSteps; ++i)
    {
        sweepA = v[edgesI[sEdge][0]] + (HReal)i*particleDiameter*normal;
        sweepB = v[edgesI[sEdge][1]] + (HReal)i*particleDiameter*normal;
        success = LineLineIntersect(v[edgesI[edge1][0]], v[edgesI[edge1][1]], sweepA, sweepB, o1, o2, m1, m2);
        i1 = o1;
        if(success == false)
        {
            std::cout << "Intersection 1 failed" << std::endl;
        }
        success = LineLineIntersect(v[edgesI[edge2][0]], v[edgesI[edge2][1]], sweepA, sweepB, o1, o2, m1, m2);
        i2 = o1;
        if(success == false)
        {
            std::cout << "Intersection 1 failed" << std::endl;
        }
        Vec3r s = i1-i2;
        int step = std::floor(s.length()/particleDiameter);
        Vec3r ps = s/((HReal)step);
        for(int j=1; j<step; ++j)
        {
            Vec3r p = i2 + (HReal)j*ps;
            samples.push_back(p);
        }
    }
    return success;
}

bool AkinciMeshSampling(const TriMesh& mesh, const HReal& particleDiameter, std::vector<Vec3r>& samples)
{
    bool success = true, tmpSuccess = false;

    //Sample vertices
    for(size_t i=0; i<mesh.vertices.size(); ++i)
        samples.push_back(mesh.vertices[i]);

    //Sample edges
    std::vector< Vec3r > tmpSample;
    std::vector< std::pair<int, int> > edges;
    mesh.getEdges(mesh.triangles, edges);
    for(size_t i=0; i<edges.size(); ++i)
    {
        tmpSuccess = AkinciEdgeSampling(mesh.vertices[edges[i].first], mesh.vertices[edges[i].second], particleDiameter, tmpSample);
        for(size_t j=0; j<tmpSample.size(); ++j)
            samples.push_back(tmpSample[j]);
        success = success && tmpSuccess;
    }

    //Sample triangles
    for(size_t i=0; i<mesh.triangles.size(); ++i)
    {
        tmpSuccess = AkinciTriangleSampling(mesh.vertices[mesh.triangles[i][0]], mesh.vertices[mesh.triangles[i][1]], mesh.vertices[mesh.triangles[i][2]], particleDiameter, tmpSample);
        for(size_t j=0; j<tmpSample.size(); ++j)
            samples.push_back(tmpSample[j]);
        success = success && tmpSuccess;
    }
    return success;
}

std::vector<Vec3r > getPyramidSampling(const Vec3r& _center, HReal base, HReal height, HReal spacing)
{
    std::vector< Vec3r > tmp;
    std::vector< Vec3r > result;
    Vec3r v1,v2,v3, center(_center[0], _center[1], _center[2]);

    //Base
    v1 = center + Vec3r(-base/2.0,0, -base/2.0);
    v2 = center + Vec3r(base/2.0, 0,-base/2.0);
    v3 = center + Vec3r(-base/2.0, 0,base/2.0);
    AkinciFullTriangleSampling(v1, v2, v3, spacing, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(Vec3r(tmp[i][0], tmp[i][1], tmp[i][2]));

    v1 = center + Vec3r(base/2.0,0, base/2.0);
    v2 = center + Vec3r(base/2.0, 0,-base/2.0);
    v3 = center + Vec3r(-base/2.0, 0,base/2.0);
    AkinciFullTriangleSampling(v1, v2, v3, spacing, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(Vec3r(tmp[i][0], tmp[i][1], tmp[i][2]));

    //Faces
    v1 = center + Vec3r(base/2.0,0, base/2.0);
    v2 = center + Vec3r(base/2.0, 0,-base/2.0);
    v3 = center + Vec3r(0, height,0);
    AkinciFullTriangleSampling(v1, v2, v3, spacing, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(Vec3r(tmp[i][0], tmp[i][1], tmp[i][2]));

    v1 = center + Vec3r(base/2.0,0, base/2.0);
    v2 = center + Vec3r(-base/2.0, 0,base/2.0);
    v3 = center + Vec3r(0, height,0);
    AkinciFullTriangleSampling(v1, v2, v3, spacing, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(Vec3r(tmp[i][0], tmp[i][1], tmp[i][2]));

    v1 = center + Vec3r(-base/2.0,0, -base/2.0);
    v2 = center + Vec3r(base/2.0, 0,-base/2.0);
    v3 = center + Vec3r(0, height,0);
    AkinciFullTriangleSampling(v1, v2, v3, spacing, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(Vec3r(tmp[i][0], tmp[i][1], tmp[i][2]));

    v1 = center + Vec3r(-base/2.0,0, -base/2.0);
    v2 = center + Vec3r(-base/2.0, 0,base/2.0);
    v3 = center + Vec3r(0, height,0);

    AkinciFullTriangleSampling(v1, v2, v3, spacing, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(Vec3r(tmp[i][0], tmp[i][1], tmp[i][2]));

    return result;
}

std::vector< Vec3r > getSphereSampling(const Vec3r& center, HReal radius, HReal spacingX, HReal spacingY)
{
    HReal theta=0.0;
    HReal phi=0.0;
    HReal l_theta = 2.0*M_PI*radius;
    HReal l_phi = M_PI*radius;
    int thetaStep = std::floor(l_theta/spacingX);
    int phiStep = std::floor(l_phi/spacingY);
    std::vector< Vec3r > result;

    for(int i=0; i<thetaStep; ++i)
    {
        theta += (2.0*M_PI/(HReal)thetaStep);
        phi=0.0;
        for(int j=0; j<phiStep; ++j)
        {
            phi += (M_PI/(HReal)phiStep);
            result.push_back( center + Vec3r(radius*cos(theta)*sin(phi), radius*sin(theta)*sin(phi), radius*cos(phi)) );
        }
    }
    return result;
}

std::vector<Vec3r > getEllipsoidSampling(const Vec3r& center, HReal axis_1, HReal axis_2, HReal axis_3, HReal spacingX, HReal spacingY)
{
    HReal theta=0.0;
    HReal phi=0.0;
    HReal l_theta = 2.0*M_PI*axis_1;
    HReal l_phi = M_PI*axis_2;
    int thetaStep = std::floor(l_theta/spacingX);
    int phiStep = std::floor(l_phi/spacingY);
    std::vector< Vec3r > result;

    for(int i=0; i<thetaStep; ++i)
    {
        theta += (2.0*M_PI/(HReal)thetaStep);
        phi=0.0;
        for(int j=0; j<phiStep; ++j)
        {
            phi += (M_PI/(HReal)phiStep);
            result.push_back( center + Vec3r(axis_1*cos(theta)*sin(phi), axis_2*sin(theta)*sin(phi), axis_3*cos(phi)) );
        }
    }
    return result;

}

std::vector<Vec3r > getCapsuleSampling(const Vec3r& center, HReal radius, HReal height, HReal spacingX, HReal spacingY)
{
    HReal theta=0.0;
    HReal phi=0.0;
    HReal l_theta = 2.0*M_PI*radius;
    HReal l_phi = (M_PI/2.0)*radius;
    int thetaStep = std::floor(l_theta/spacingX);
    int phiStep = std::floor(l_phi/spacingY);
    std::vector< Vec3r > result;
    Vec3r c1(center[0], center[2], height);
    Vec3r c2(center[0], center[2], 0);

    for(int i=0; i<thetaStep; ++i)
    {
        theta += (2.0*M_PI/(HReal)thetaStep);
        phi=0.0;
        for(int j=0; j<phiStep; ++j)
        {
            phi += (M_PI/(HReal)(2.0*phiStep));
            result.push_back( c1 + Vec3r(radius*cos(theta)*sin(phi), radius*sin(theta)*sin(phi), radius*cos(phi)) );
        }
    }

    for(int i=0; i<thetaStep; ++i)
    {
        theta += (2.0*M_PI/(HReal)thetaStep);
        phi=M_PI/2.0-(M_PI/(2.0*phiStep));
        for(int j=0; j<phiStep; ++j)
        {
            phi += (M_PI/(HReal)(2.0*(phiStep)));
            result.push_back( c2 + Vec3r(radius*cos(theta)*sin(phi), radius*sin(theta)*sin(phi), radius*cos(phi)) );
        }
    }

    std::vector< Vec3r > tmp = getCylinderSampling(center, height, radius, spacingX, spacingY);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(tmp[i]);

    return result;
}

std::vector< Vec3r > getHemiSphereSampling(const Vec3r& center, HReal radius, HReal spacingX, HReal spacingY)
{
    HReal theta=0.0;
    HReal phi=0.0;
    HReal l_theta = 2.0*M_PI*radius;
    HReal l_phi = (M_PI/2.0)*radius;
    int thetaStep = std::floor(l_theta/spacingX);
    int phiStep = std::floor(l_phi/spacingY);
    std::vector< Vec3r > result;

    for(int i=0; i<thetaStep; ++i)
    {
        theta += (2.0*M_PI/(HReal)thetaStep);
        phi=0.0;
        for(int j=0; j<phiStep; ++j)
        {
            phi += (M_PI/(HReal)(2.0*phiStep));
            result.push_back( center + Vec3r(radius*cos(theta)*sin(phi), radius*sin(theta)*sin(phi), radius*cos(phi)) );
        }
    }
    return result;
}

std::vector< Vec3r > getTorusSampling(const Vec3r& center, HReal tubeRadius, HReal innerRadius, HReal spacingX, HReal spacingY)
{
    HReal u=0.0;
    HReal v=0.0;
    HReal l_u = 2.0*M_PI*innerRadius;
    HReal l_v = 2.0*M_PI*tubeRadius;
    int uStep = std::floor(l_u/spacingX);
    int vStep = std::floor(l_v/spacingY);
    std::vector< Vec3r > result;

    for(int i=0; i<uStep; ++i)
    {
        u += (2.0*M_PI/(HReal)uStep);
        v =0.0;
        for(int j=0; j<vStep; ++j)
        {
            v += (2.0*M_PI/(HReal)vStep);
            result.push_back( center + Vec3r( (innerRadius+tubeRadius*cos(v))*cos(u), (innerRadius+tubeRadius*cos(v))*sin(u), tubeRadius*sin(v)) );
        }
    }
    return result;
}

std::vector< Vec3r > getConeSampling(const Vec3r& center, HReal height, HReal stopHeight, HReal baseRadius, HReal spacingX, HReal spacingY)
{
    HReal theta=0.0;
    HReal u=0.0;
    HReal l_theta = 2.0*M_PI*baseRadius;
    int thetaStep = std::floor(l_theta/spacingX);
    int uStep = std::floor(stopHeight/spacingY);
    std::vector< Vec3r > result;

    for(int i=0; i<thetaStep; ++i)
    {
        theta += (2.0*M_PI/(HReal)thetaStep);
        u =0.0;
        for(int j=0; j<uStep; ++j)
        {
            u += (stopHeight/(HReal)uStep);
            result.push_back( center + Vec3r( ((height-u)/height)*baseRadius*cos(theta), ((height-u)/height)*baseRadius*sin(theta), u));
        }
    }
    return result;
}


std::vector< Vec3r > getCylinderSampling(const Vec3r& center, HReal height, HReal baseRadius, HReal spacingX, HReal spacingY)
{
    HReal theta=0.0;
    HReal u=0.0;
    HReal l_theta = 2.0*M_PI*baseRadius;
    int thetaStep = std::floor(l_theta/spacingX);
    int uStep = std::floor(height/spacingY);
    std::vector< Vec3r > result;

    for(int i=0; i<thetaStep; ++i)
    {
        theta += (2.0*M_PI/(HReal)thetaStep);
        u =0.0;
        for(int j=0; j<uStep; ++j)
        {
            u += (height/(HReal)uStep);
            result.push_back( center + Vec3r( baseRadius*cos(theta), baseRadius*sin(theta), u));
        }
    }
    return result;
}

std::vector<Vec3r> getCubeSampling(const Vec3r& offset, const Vec3r& scale, HReal spacing)
{
    std::vector<Vec3r> result;
    int widthSize = floor(scale[0]/spacing);
    int heightSize = floor(scale[1]/spacing);
    int depthSize = floor(scale[2]/spacing);

    for(int i=0; i<widthSize; ++i)
    {
        for(int j=0; j<heightSize; ++j)
        {
            for(int k=0; k<depthSize; ++k)
            {
                result.push_back(offset + Vec3r(i*spacing,j*spacing,k*spacing));
            }
        }
    }

    return result;
}

std::vector<Vec3r> getBoxSampling(const Vec3r& offset, const Vec3r& scale, HReal spacing)
{
    int epsilon = 0;
    int widthSize = std::floor(scale[0]/spacing);
    int heightSize = std::floor(scale[1]/spacing);
    int depthSize = std::floor(scale[2]/spacing);
    std::vector<Vec3r> result;

    //ZX plane - bottom
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            result.push_back(Vec3r(offset[0]+i*spacing, offset[1], offset[2]+j*spacing));
        }
    }

    //ZX plane - top
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            result.push_back(Vec3r(offset[0]+i*spacing, offset[1]+scale[1], offset[2]+j*spacing));
        }
    }

    //XY plane - back
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= heightSize+epsilon; ++j)
        {
            result.push_back(Vec3r(offset[0]+i*spacing, offset[1]+j*spacing, offset[2]));
        }
    }

    //XY plane - front
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= heightSize-epsilon; ++j)
        {
            result.push_back(Vec3r(offset[0]+i*spacing, offset[1]+j*spacing, offset[2]+scale[2]));
        }
    }

    //YZ plane - left
    for(int i = -epsilon; i <= heightSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            result.push_back(Vec3r(offset[0], offset[1]+i*spacing, offset[2]+j*spacing));
        }
    }

    //YZ plane - right
    for(int i = -epsilon; i <= heightSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            result.push_back(Vec3r(offset[0]+scale[0], offset[1]+i*spacing, offset[2]+j*spacing));
        }
    }
    return result;
}

std::vector<Vec3r> getBallSampling(const Vec3r& center, HReal radius, HReal spacing)
{
    std::vector<Vec3r> result;
    Vec3r scale(2.0*radius, 2.0*radius, 2.0*radius);
    Vec3r offset = center - Vec3r(radius, radius, radius);
    GridUtility grid(offset, scale, spacing);

    for(int i=0; i<grid.size(); ++i)
    {
        Vec3r x = grid.gridToWorld(i);
        x+=grid.spacing()/2.0;
        HReal l2 = (center-x).lengthSquared();
        if(l2<=(radius*radius))
        {
            result.push_back(x);
        }
    }
    return result;
}

std::vector<Vec3r > getDiskSampling(const Vec3r& center, HReal radius, HReal spacing)
{
    HReal theta=0.0;
    HReal u=0.0;
    HReal l_theta = 2.0*M_PI*radius;
    int thetaStep = std::floor(l_theta/spacing);
    int uStep = std::floor(radius/spacing);
    std::vector< Vec3r > result;

    for(int i=0; i<thetaStep; ++i)
    {
        theta += (2.0*M_PI/(HReal)thetaStep);
        u =0.0;
        for(int j=0; j<uStep; ++j)
        {
            u += (radius/(HReal)uStep);
            result.push_back( center + Vec3r( u*radius*cos(theta), u*radius*sin(theta), 0));
        }
    }
    return result;
}

}
