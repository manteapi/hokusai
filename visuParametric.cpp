#ifndef VISUPARAMETRIC_CPP
#define VISUPARAMETRIC_CPP

#include "./../include/visuParametric.hpp"
#include "./../include/sampler.hpp"
#include "./../include/collision.hpp"
#include <cmath>
#include <iostream>
#include "./../include/io.hpp"

using namespace std;
using namespace glm;

std::vector<Vec3<double> > VisuParametric::getPyramidSampling(const Vec3<double>& _center, double base, double height, double spacing)
{
    std::vector< Vec3<float> > tmp;
    std::vector< Vec3<double> > result;
    Vec3<float> v1,v2,v3,v4,v5, center(_center[0], _center[1], _center[2]);

    //Base
    v1 = center + Vec3<float>(-base/2.0,0, -base/2.0);
    v2 = center + Vec3<float>(base/2.0, 0,-base/2.0);
    v3 = center + Vec3<float>(-base/2.0, 0,base/2.0);
    AkinciFullTriangleSampling(v1, v2, v3, spacing, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(Vec3<double>(tmp[i][0], tmp[i][1], tmp[i][2]));

    v1 = center + Vec3<float>(base/2.0,0, base/2.0);
    v2 = center + Vec3<float>(base/2.0, 0,-base/2.0);
    v3 = center + Vec3<float>(-base/2.0, 0,base/2.0);
    AkinciFullTriangleSampling(v1, v2, v3, spacing, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(Vec3<double>(tmp[i][0], tmp[i][1], tmp[i][2]));

    //Faces
    v1 = center + Vec3<float>(base/2.0,0, base/2.0);
    v2 = center + Vec3<float>(base/2.0, 0,-base/2.0);
    v3 = center + Vec3<float>(0, height,0);
    AkinciFullTriangleSampling(v1, v2, v3, spacing, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(Vec3<double>(tmp[i][0], tmp[i][1], tmp[i][2]));

    v1 = center + Vec3<float>(base/2.0,0, base/2.0);
    v2 = center + Vec3<float>(-base/2.0, 0,base/2.0);
    v3 = center + Vec3<float>(0, height,0);
    AkinciFullTriangleSampling(v1, v2, v3, spacing, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(Vec3<double>(tmp[i][0], tmp[i][1], tmp[i][2]));

    v1 = center + Vec3<float>(-base/2.0,0, -base/2.0);
    v2 = center + Vec3<float>(base/2.0, 0,-base/2.0);
    v3 = center + Vec3<float>(0, height,0);
    AkinciFullTriangleSampling(v1, v2, v3, spacing, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(Vec3<double>(tmp[i][0], tmp[i][1], tmp[i][2]));

    v1 = center + Vec3<float>(-base/2.0,0, -base/2.0);
    v2 = center + Vec3<float>(-base/2.0, 0,base/2.0);
    v3 = center + Vec3<float>(0, height,0);

    AkinciFullTriangleSampling(v1, v2, v3, spacing, tmp);
    for(size_t i=0; i<tmp.size(); ++i)
        result.push_back(Vec3<double>(tmp[i][0], tmp[i][1], tmp[i][2]));

    return result;
}

std::vector< Vec3<double> > VisuParametric::getSphereSampling(const Vec3<double>& center, double radius, double spacingX, double spacingY)
{
    double theta=0.0;
    double phi=0.0;
    double l_theta = 2.0*M_PI*radius;
    double l_phi = M_PI*radius;
    int thetaStep = std::floor(l_theta/spacingX);
    int phiStep = std::floor(l_phi/spacingY);
    std::vector< Vec3<double> > result;

    for(int i=0; i<thetaStep; ++i)
    {
        theta += (2.0*M_PI/(double)thetaStep);
        phi=0.0;
        for(int j=0; j<phiStep; ++j)
        {
            phi += (M_PI/(double)phiStep);
            result.push_back( center + Vec3<double>(radius*cos(theta)*sin(phi), radius*sin(theta)*sin(phi), radius*cos(phi)) );
        }
    }
    return result;
}

    std::vector<Vec3<double> > VisuParametric::getEllipsoidSampling(const Vec3<double>& center, double axis_1, double axis_2, double axis_3, double spacingX, double spacingY)
    {
        double theta=0.0;
        double phi=0.0;
        double l_theta = 2.0*M_PI*axis_1;
        double l_phi = M_PI*axis_2;
        int thetaStep = std::floor(l_theta/spacingX);
        int phiStep = std::floor(l_phi/spacingY);
        std::vector< Vec3<double> > result;

        for(int i=0; i<thetaStep; ++i)
        {
            theta += (2.0*M_PI/(double)thetaStep);
            phi=0.0;
            for(int j=0; j<phiStep; ++j)
            {
                phi += (M_PI/(double)phiStep);
                result.push_back( center + Vec3<double>(axis_1*cos(theta)*sin(phi), axis_2*sin(theta)*sin(phi), axis_3*cos(phi)) );
            }
        }
        return result;

    }

    std::vector<Vec3<double> > VisuParametric::getCapsuleSampling(const Vec3<double>& center, double radius, double height, double spacingX, double spacingY)
    {
        double theta=0.0;
        double phi=0.0;
        double l_theta = 2.0*M_PI*radius;
        double l_phi = (M_PI/2.0)*radius;
        int thetaStep = std::floor(l_theta/spacingX);
        int phiStep = std::floor(l_phi/spacingY);
        std::vector< Vec3<double> > result;
        Vec3<double> c1(center[0], center[2], height);
        Vec3<double> c2(center[0], center[2], 0);

                for(int i=0; i<thetaStep; ++i)
        {
            theta += (2.0*M_PI/(double)thetaStep);
            phi=0.0;
            for(int j=0; j<phiStep; ++j)
            {
                phi += (M_PI/(double)(2.0*phiStep));
                result.push_back( c1 + Vec3<double>(radius*cos(theta)*sin(phi), radius*sin(theta)*sin(phi), radius*cos(phi)) );
            }
        }

        for(int i=0; i<thetaStep; ++i)
        {
            theta += (2.0*M_PI/(double)thetaStep);
            phi=M_PI/2.0-(M_PI/(2.0*phiStep));
            for(int j=0; j<phiStep; ++j)
            {
                phi += (M_PI/(double)(2.0*(phiStep)));
                result.push_back( c2 + Vec3<double>(radius*cos(theta)*sin(phi), radius*sin(theta)*sin(phi), radius*cos(phi)) );
            }
        }

        std::vector< Vec3<double> > tmp = getCylinderSampling(center, height, radius, spacingX, spacingY);
        for(size_t i=0; i<tmp.size(); ++i)
            result.push_back(tmp[i]);

        return result;
    }

std::vector< Vec3<double> > VisuParametric::getHemiSphereSampling(const Vec3<double>& center, double radius, double spacingX, double spacingY)
{
    double theta=0.0;
    double phi=0.0;
    double l_theta = 2.0*M_PI*radius;
    double l_phi = (M_PI/2.0)*radius;
    int thetaStep = std::floor(l_theta/spacingX);
    int phiStep = std::floor(l_phi/spacingY);
    std::vector< Vec3<double> > result;

    for(int i=0; i<thetaStep; ++i)
    {
        theta += (2.0*M_PI/(double)thetaStep);
        phi=0.0;
        for(int j=0; j<phiStep; ++j)
        {
            phi += (M_PI/(double)(2.0*phiStep));
            result.push_back( center + Vec3<double>(radius*cos(theta)*sin(phi), radius*sin(theta)*sin(phi), radius*cos(phi)) );
        }
    }
    return result;
}

std::vector< Vec3<double> > VisuParametric::getTorusSampling(const Vec3<double>& center, double tubeRadius, double innerRadius, double spacingX, double spacingY)
{
    double u=0.0;
    double v=0.0;
    double l_u = 2.0*M_PI*innerRadius;
    double l_v = 2.0*M_PI*tubeRadius;
    int uStep = std::floor(l_u/spacingX);
    int vStep = std::floor(l_v/spacingY);
    std::vector< Vec3<double> > result;

    for(int i=0; i<uStep; ++i)
    {
        u += (2.0*M_PI/(double)uStep);
        v =0.0;
        for(int j=0; j<vStep; ++j)
        {
            v += (2.0*M_PI/(double)vStep);
            result.push_back( center + Vec3<double>( (innerRadius+tubeRadius*cos(v))*cos(u), (innerRadius+tubeRadius*cos(v))*sin(u), tubeRadius*sin(v)) );
        }
    }
    return result;
}

std::vector< Vec3<double> > VisuParametric::getConeSampling(const Vec3<double>& center, double height, double stopHeight, double baseRadius, double spacingX, double spacingY)
{
    double theta=0.0;
    double u=0.0;
    double l_theta = 2.0*M_PI*baseRadius;
    int thetaStep = std::floor(l_theta/spacingX);
    int uStep = std::floor(stopHeight/spacingY);
    std::vector< Vec3<double> > result;

    for(int i=0; i<thetaStep; ++i)
    {
        theta += (2.0*M_PI/(double)thetaStep);
        u =0.0;
        for(int j=0; j<uStep; ++j)
        {
            u += (stopHeight/(double)uStep);
            result.push_back( center + Vec3<double>( ((height-u)/height)*baseRadius*cos(theta), ((height-u)/height)*baseRadius*sin(theta), u));
        }
    }
    return result;
}


std::vector< Vec3<double> > VisuParametric::getCylinderSampling(const Vec3<double>& center, double height, double baseRadius, double spacingX, double spacingY)
{
    double theta=0.0;
    double u=0.0;
    double l_theta = 2.0*M_PI*baseRadius;
    int thetaStep = std::floor(l_theta/spacingX);
    int uStep = std::floor(height/spacingY);
    std::vector< Vec3<double> > result;

    for(int i=0; i<thetaStep; ++i)
    {
        theta += (2.0*M_PI/(double)thetaStep);
        u =0.0;
        for(int j=0; j<uStep; ++j)
        {
            u += (height/(double)uStep);
            result.push_back( center + Vec3<double>( baseRadius*cos(theta), baseRadius*sin(theta), u));
        }
    }
    return result;
}

std::vector<Vec3<double> > VisuParametric::getDiskSampling(const Vec3<double>& center, double radius, double spacing)
{
    double theta=0.0;
    double u=0.0;
    double l_theta = 2.0*M_PI*radius;
    int thetaStep = std::floor(l_theta/spacing);
    int uStep = std::floor(radius/spacing);
    std::vector< Vec3<double> > result;

    for(int i=0; i<thetaStep; ++i)
    {
        theta += (2.0*M_PI/(double)thetaStep);
        u =0.0;
        for(int j=0; j<uStep; ++j)
        {
            u += (radius/(double)uStep);
            result.push_back( center + Vec3<double>( u*radius*cos(theta), u*radius*sin(theta), 0));
        }
    }
    return result;
}

VisuParametric::VisuParametric()
{
    picked = 0;
    model = mat4(1.0);
    radius = 0.0;
    center = vec3(0,0,0);
    std::cout << "out" << std::endl;

    Vec3<double> center(0,0,0);
    double radius = 2.0, spacingX = 0.1, spacingY = 0.1;
    double axis_1 = 1.0, axis_2 = 1.0, axis_3 = 1.0;
    //std::vector< Vec3<double> > tmpPoints = getSphereSampling(center, radius, spacingX, spacingY);
    //std::vector< Vec3<double> > tmpPoints = getHemiSphereSampling(center, radius, spacingX, spacingY);
    //double tubeRadius = 0.5, innerRadius = 1.0;
    //std::vector< Vec3<double> > tmpPoints = getTorusSampling(center, tubeRadius, innerRadius, spacingX, spacingY);

    double baseRadius = 0.5; double height = 1.0;
//    std::vector< Vec3<double> > tmpPoints = getHoofSampling(center, radius, height, spacingX);
//    std::vector< Vec3<double> > tmpPoints = getCapsuleSampling(center, radius, height, spacingX, spacingY);
//    std::vector< Vec3<double> > tmpPoints = getConeSampling(center, height, height/2.0, baseRadius, spacingX, spacingY);
//    std::vector< Vec3<double> > tmpPoints = getConeSampling(center, height, baseRadius, spacingX, spacingY);
    std::vector< Vec3<double> > tmpPoints = getPyramidSampling(center, radius, height, spacingX);
    //std::vector< Vec3<double> > tmpPoints = getDiskSampling(center, radius, spacingX);

    for(size_t i=0; i<tmpPoints.size(); ++i)
    {
        vertices.push_back(getVec(tmpPoints[i]));
        colors.push_back(vec4(1,0,0,1));
    }
    std::cout << vertices.size() << std::endl;
}

VisuParametric::~VisuParametric()
{
    glDeleteBuffers(1, &positionBuffer);
    glDeleteBuffers(1, &colorBuffer);
}

void VisuParametric::initShader( int shaderId )
{
    modelLoc = glGetUniformLocation(shaderId, "modelMat");
    positionLoc = glGetAttribLocation(shaderId, "position");
    colorLoc = glGetAttribLocation(shaderId, "vertex_color");
}

void VisuParametric::createvbo()
{
    //VisuParametric line data
    glGenBuffers(1, &positionBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, positionBuffer);
    glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(glm::vec3), vertices.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &colorBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, colorBuffer);
    glBufferData(GL_ARRAY_BUFFER, colors.size()*sizeof(glm::vec4), colors.data(), GL_STATIC_DRAW);
}

void VisuParametric::reloadvbo()
{
    //VisuParametric line data
    glBindBuffer(GL_ARRAY_BUFFER, positionBuffer);
    glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(glm::vec3), vertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, colorBuffer);
    glBufferData(GL_ARRAY_BUFFER, colors.size()*sizeof(glm::vec4), colors.data(), GL_STATIC_DRAW);
}

void VisuParametric::drawVisuParametric()
{
    //Send uniform
    this->model = mat4(1.0f);
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

    //Active/Link vertices
    glEnableVertexAttribArray(positionLoc);
    glBindBuffer(GL_ARRAY_BUFFER,positionBuffer);
    glVertexAttribPointer(positionLoc,3,GL_FLOAT,GL_FALSE,0,(void *)0);

    glEnableVertexAttribArray(colorLoc);
    glBindBuffer(GL_ARRAY_BUFFER,colorBuffer);
    glVertexAttribPointer(colorLoc,4,GL_FLOAT,GL_FALSE,0,(void *)0);

    glDrawArrays(GL_POINTS, 0, vertices.size());

    glDisableVertexAttribArray(positionLoc);
    glDisableVertexAttribArray(colorLoc);
}

void VisuParametric::drawShader()
{
    drawVisuParametric();
}

void VisuParametric::drawFix(){}

float VisuParametric::rayPicking(vec3& source, vec3& ray)
{}

void VisuParametric::keyPressedEvent(sf::Event& e )
{}
#endif
