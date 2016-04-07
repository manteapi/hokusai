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

#ifndef HOKUSAI_SYSTEM_HPP
#define HOKUSAI_SYSTEM_HPP

#include <iostream>
#include <fstream>
#include <ctime>

#include <aljabr/AljabrCore>
#include "utils.hpp"
#include "common.hpp"
#include "particle.hpp"
#include "gridUtility.hpp"
#include "triMesh.hpp"
#include "sampler.hpp"
#include "particleSource.hpp"

namespace hokusai
{

class System
{

public:

    System();
    System(int resolution);
    ~System();

public :
    int countTime;
    int countExport;
    int particleNumber;
    int boundaryNumber;

    HReal particlePerCell;
    HReal volume;
    HReal restDensity;
    HReal mean_density;
    HReal density_fluctuation;
    HReal real_volume;
    HReal mass;
    HReal h; // Smoothing radius
    HReal fcohesion; //Fluid cohesion
    HReal badhesion; //Boundary adhesion
    HReal sigma; //Boundary friction
    HReal cs;// Sound speed
    HReal alpha; // Viscosity
    HReal boundaryH;
    HReal dt;
    HReal time;
    HReal rho_avg_l;
    HReal maxEta;

    Vec3r gravity;

    AkinciKernel a_kernel;
    MonaghanKernel p_kernel;
    BoundaryKernel b_kernel;

    vector<Particle> particles;
    vector<Boundary> boundaries;

    GridUtility gridInfo;
    vector< vector<int> > boundaryGrid;
    vector< vector<int> > fluidGrid;

    vector< ParticleSource > p_sources;

public :
    void getNearestNeighbor(vector< int >& neighbors, const vector<vector<int> > &grid, const Vec3r& x);
    void getNearestNeighbor(const int i, const HReal radius);

    //Simulation Loop
    void prepareGrid();
    void computeSurfaceParticle();
    void predictAdvection();
    void predictRho(int i);
    void initializePressure(int i);
    void computeNormal(int i);
    bool isSurfaceParticle(int i, HReal treshold);
    vector<Particle> getSurfaceParticle();
    void computeRho(int i);
    void computeAdvectionForces(int i);
    void predictVelocity(int i);
    void computeDii(int i);
    void computeDii_Fluid(int i);
    void computeDii_Boundary(int i);
    void computeAii(int i);
    void pressureSolve();
    void computeSumDijPj(int i);
    
    void computeViscosityForces(int i, int j)
    {
        Particle& pi=particles[i];
        Particle& pj=particles[j];
        Vec3r r = pi.x - pj.x;
        Vec3r vij = pi.v - pj.v;
        HReal dotVijRij = Vec3r::dotProduct(vij,r);
        if(dotVijRij < 0)
        {
            HReal kij = 2.0*restDensity/(pi.rho+pj.rho);
            HReal epsilon=0.01;
            Vec3r gradient(0.0);
            p_kernel.monaghanGradient(r, gradient);
            HReal Pij = -kij*(2.0*alpha*h*cs/(pi.rho+pj.rho)) * ( dotVijRij / (r.lengthSquared() + epsilon*h*h) );
            pi.f_adv += -kij*mass*mass*Pij*gradient;
        }
    }

    void computeBoundaryFrictionForces(int i, int j)
    {
        Particle& pi=particles[i];
            Boundary& bj=boundaries[j];
            Vec3r vij = pi.v;//-pj.v;
            Vec3r xij= pi.x - bj.x;
            HReal dotVijRij = Vec3r::dotProduct(vij,xij);
            if(dotVijRij<0)
            {
                Vec3r gradient(0.0);
                HReal epsilon=0.01;
                HReal nu = (sigma*h*cs)/(2.0*pi.rho);
                HReal Pij = -nu * ( std::min(dotVijRij,0.0) / (xij.lengthSquared() + epsilon*h*h) );
                p_kernel.monaghanGradient(xij, gradient);
                pi.f_adv += -mass*bj.psi*Pij*gradient;
            }
    }

    void computeSurfaceTensionForces(int i, int j)
    {
        if(i!=j)
        {
            Particle& pi=particles[i];
            Particle& pj=particles[j];
            if(pi.isSurface==true || pj.isSurface==true)
            {
                Vec3r r = pi.x - pj.x;
                HReal kij = 2.0*restDensity/(pi.rho+pj.rho);
                HReal l = r.length();
                Vec3r cohesionForce = -(fcohesion*mass*mass*a_kernel.cohesionValue(l)/l) * r;
                Vec3r nij = pi.n-pj.n;
                Vec3r curvatureForce = -fcohesion*mass*nij;
                pi.f_adv += kij*(cohesionForce+curvatureForce);
            }
        }
    }
    
    void computeBoundaryAdhesionForces(int i, int j)
    {
        Particle& pi=particles[i];
            Boundary& bj=boundaries[j];
            Vec3r xij= pi.x - bj.x;
            HReal l = xij.length();
            pi.f_adv += -(badhesion*mass*boundaries[j].psi*a_kernel.adhesionValue(l)/l)*xij;
    }

    Vec3r computeDij(int i, int j)
    {
        Particle& pi=particles[i];
        Particle& pj=particles[j];
        Vec3r gradient(0.0);
        p_kernel.monaghanGradient(pi.x-pj.x, gradient);
        Vec3r d=-(dt*dt*mass)/pow(pj.rho,2)*gradient;
        return d;
    }

    void computePressure(int i);
    void computePressureForce(int i);

    void computeFluidPressureForce(int i, int j)
    {
        Vec3r gradient(0.0);
        Particle& pi=particles[i];
        Particle& pj=particles[j];
        p_kernel.monaghanGradient(pi.x-pj.x, gradient);
        if( i!=j )
        {
            pi.f_p += -mass*mass*( pi.p/pow(pi.rho,2) + pj.p/pow(pj.rho,2) ) * gradient;
        }
    }

    void computeBoundaryPressureForce(int i, int j)
    {
        Vec3r gradient(0.0);
        Particle& pi=particles[i];
        Boundary& bj=boundaries[j];
            p_kernel.monaghanGradient(pi.x-bj.x, gradient);
            pi.f_p += -mass*bj.psi*( pi.p/pow(pi.rho,2) ) * gradient;
    }

    void computeError();
    void integration();
    void simulate();
    void applySources();
    void applySinks();

    void addBoundaryParticle(const Vec3r& x, const Vec3r& v = Vec3r(0,0,0));
    void addFluidParticle(const Vec3r& x, const Vec3r& v = Vec3r(0,0,0));

    //Initialize a dam break scenario
    const Vec3r& getGravity();
    void setGravity(const Vec3r& _gravity);
    void setParameters(int _number, HReal _volume=1.0);
    void createParticleVolume(Vec3r& pos, HReal width, HReal height, HReal depth, HReal spacing, int particleMax);

    void translateParticles(const Vec3r& t);
    void translateBoundaries(const Vec3r& t);

    void addParticleBox(HReal width, HReal height, HReal depth, HReal spacing);
    void addParticleBox(const Vec3r& offset, const Vec3r& dimension);
    void addParticleSphere(const Vec3r& centre, const HReal radius);

    void addParticleSource(const ParticleSource& s);

    //Boundary sampling
    void addBoundaryMesh(const char* filename);
    void addBoundaryBox(const Vec3r& offset, const Vec3r& scale);
    void addBoundarySphere(const Vec3r& offset, const HReal& radius);
    void addBoundaryHemiSphere(const Vec3r& offset, const HReal& radius);
    void addBoundaryDisk(const Vec3r& offset, const HReal& radius);

    void debugFluid();
    void debugIteration(int l);
    void mortonSort();
    void init();

    //Data
    void computeBoundaryVolume();
    void computeMeanDensity();
    void computeStats();
    void computeDensityFluctuation();
    void computeVolume();
    void computeCellChange();
    void computeScalarField();
    void computeSurface();

    //Getter
    vector< Vec3r > getPosition(){ vector<Vec3r > pos; for(int i=0; i<particleNumber; ++i){pos.push_back(particles[i].x);} return pos;}
    vector< Vec3r > getVelocity(){ vector<Vec3r > vel; for(int i=0; i<particleNumber; ++i){vel.push_back(particles[i].v);} return vel;}
    vector< Vec3r > getNormal(){ vector<Vec3r > normal; for(int i=0; i<particleNumber; ++i){normal.push_back(particles[i].n);} return normal;}
    vector< HReal > getDensity(){ vector<HReal> density; for(int i=0; i<particleNumber; ++i){density.push_back(particles[i].rho);} return density;}
    vector< HReal > getMass(){ vector<HReal> o_mass; for(int i=0; i<particleNumber; ++i){o_mass.push_back(mass);} return o_mass;}

    void write(const char* filename, vector< Vec3r > data);
    void write(const char* filename, vector<HReal> data);
    void exportState(const char* baseName);
    HReal getTime(){return time;}
    HReal & getSmoothingRadiusValue(){ return h; }
    const HReal & getSmoothingRadius() const { return h; }
    HReal & getTimeStepValue(){ return dt; }
    HReal & getMassValue(){ return mass; }
    HReal & getMeanDensityValue(){ return mean_density;}
    HReal & getDensityFluctuationValue(){ return density_fluctuation;}
    HReal & getRealVolumeValue(){ return real_volume; }
    int & getParticleNumber(){ return particleNumber; }

    HReal & getViscosity(){return alpha;}
    const HReal & getViscosity() const {return alpha;}
    HReal & getFluidCohesion(){return fcohesion;}
    const HReal & getFluidCohesion() const {return fcohesion;}
    HReal & getBoundaryAdhesion(){return badhesion;}
    const HReal & getBoundaryAdhesion() const {return badhesion;}
    const HReal & getBoundaryFriction() const {return sigma;}
    HReal & getBoundaryFriction() {return sigma;}
    HReal & getTimeStep() {return dt;}
    const HReal& getTimeStep() const {return dt;}
};

bool pairCompare( const std::pair<int,int>& e1, const std::pair<int,int>& e2 );
}

#endif // SYSTEM_H
