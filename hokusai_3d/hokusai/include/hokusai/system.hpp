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
    int m_countTime;
    int m_countExport;
    int m_particleNumber;
    int m_boundaryNumber;

    HReal m_particlePerCell;
    HReal m_volume;
    HReal m_restDensity;
    HReal m_meanDensity;
    HReal m_densityFluctuation;
    HReal m_realVolume;
    HReal m_mass;
    HReal m_h; // Smoothing radius
    HReal m_fcohesion; //Fluid cohesion
    HReal m_badhesion; //Boundary adhesion
    HReal m_sigma; //Boundary friction
    HReal m_cs;// Sound speed
    HReal m_alpha; // Viscosity
    HReal m_boundaryH;
    HReal m_dt;
    HReal m_time;
    HReal m_rho_avg_l;
    HReal m_maxEta;

    Vec3r m_gravity;

    AkinciKernel m_aKernel;
    MonaghanKernel m_pKernel;
    BoundaryKernel m_bKernel;

    vector<Particle> m_particles;
    vector<Boundary> m_boundaries;

    GridUtility m_gridInfo;
    vector< vector<int> > m_boundaryGrid;
    vector< vector<int> > m_fluidGrid;

    vector< ParticleSource > m_pSources;

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
        Particle& pi=m_particles[i];
        Particle& pj=m_particles[j];
        Vec3r r = pi.x - pj.x;
        Vec3r vij = pi.v - pj.v;
        HReal dotVijRij = Vec3r::dotProduct(vij,r);
        if(dotVijRij < 0)
        {
            HReal kij = 2.0*m_restDensity/(pi.rho+pj.rho);
            HReal epsilon=0.01;
            Vec3r gradient(0.0);
            m_pKernel.monaghanGradient(r, gradient);
            HReal Pij = -kij*(2.0*m_alpha*m_h*m_cs/(pi.rho+pj.rho)) * ( dotVijRij / (r.lengthSquared() + epsilon*m_h*m_h) );
            pi.f_adv += -kij*m_mass*m_mass*Pij*gradient;
        }
    }

    void computeBoundaryFrictionForces(int i, int j)
    {
        Particle& pi=m_particles[i];
            Boundary& bj=m_boundaries[j];
            Vec3r vij = pi.v;//-pj.v;
            Vec3r xij= pi.x - bj.x;
            HReal dotVijRij = Vec3r::dotProduct(vij,xij);
            if(dotVijRij<0)
            {
                Vec3r gradient(0.0);
                HReal epsilon=0.01;
                HReal nu = (m_sigma*m_h*m_cs)/(2.0*pi.rho);
                HReal Pij = -nu * ( std::min(dotVijRij,0.0) / (xij.lengthSquared() + epsilon*m_h*m_h) );
                m_pKernel.monaghanGradient(xij, gradient);
                pi.f_adv += -m_mass*bj.psi*Pij*gradient;
            }
    }

    void computeSurfaceTensionForces(int i, int j)
    {
        if(i!=j)
        {
            Particle& pi=m_particles[i];
            Particle& pj=m_particles[j];
            if(pi.isSurface==true || pj.isSurface==true)
            {
                Vec3r r = pi.x - pj.x;
                HReal kij = 2.0*m_restDensity/(pi.rho+pj.rho);
                HReal l = r.length();
                Vec3r cohesionForce = -(m_fcohesion*m_mass*m_mass*m_aKernel.cohesionValue(l)/l) * r;
                Vec3r nij = pi.n-pj.n;
                Vec3r curvatureForce = -m_fcohesion*m_mass*nij;
                pi.f_adv += kij*(cohesionForce+curvatureForce);
            }
        }
    }
    
    void computeBoundaryAdhesionForces(int i, int j)
    {
        Particle& pi=m_particles[i];
            Boundary& bj=m_boundaries[j];
            Vec3r xij= pi.x - bj.x;
            HReal l = xij.length();
            pi.f_adv += -(m_badhesion*m_mass*m_boundaries[j].psi*m_aKernel.adhesionValue(l)/l)*xij;
    }

    Vec3r computeDij(int i, int j)
    {
        Particle& pi=m_particles[i];
        Particle& pj=m_particles[j];
        Vec3r gradient(0.0);
        m_pKernel.monaghanGradient(pi.x-pj.x, gradient);
        Vec3r d=-(m_dt*m_dt*m_mass)/pow(pj.rho,2)*gradient;
        return d;
    }

    void computePressure(int i);
    void computePressureForce(int i);

    void computeFluidPressureForce(int i, int j)
    {
        Vec3r gradient(0.0);
        Particle& pi=m_particles[i];
        Particle& pj=m_particles[j];
        m_pKernel.monaghanGradient(pi.x-pj.x, gradient);
        if( i!=j )
        {
            pi.f_p += -m_mass*m_mass*( pi.p/pow(pi.rho,2) + pj.p/pow(pj.rho,2) ) * gradient;
        }
    }

    void computeBoundaryPressureForce(int i, int j)
    {
        Vec3r gradient(0.0);
        Particle& pi=m_particles[i];
        Boundary& bj=m_boundaries[j];
            m_pKernel.monaghanGradient(pi.x-bj.x, gradient);
            pi.f_p += -m_mass*bj.psi*( pi.p/pow(pi.rho,2) ) * gradient;
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
    vector< Vec3r > getPosition(){ vector<Vec3r > pos; for(int i=0; i<m_particleNumber; ++i){pos.push_back(m_particles[i].x);} return pos;}
    vector< Vec3r > getVelocity(){ vector<Vec3r > vel; for(int i=0; i<m_particleNumber; ++i){vel.push_back(m_particles[i].v);} return vel;}
    vector< Vec3r > getNormal(){ vector<Vec3r > normal; for(int i=0; i<m_particleNumber; ++i){normal.push_back(m_particles[i].n);} return normal;}
    vector< HReal > getDensity(){ vector<HReal> density; for(int i=0; i<m_particleNumber; ++i){density.push_back(m_particles[i].rho);} return density;}
    vector< HReal > getMass(){ vector<HReal> o_mass; for(int i=0; i<m_particleNumber; ++i){o_mass.push_back(m_mass);} return o_mass;}

    void write(const char* filename, vector< Vec3r > data);
    void write(const char* filename, vector<HReal> data);
    void exportState(const char* baseName);
    HReal getTime(){return m_time;}
    HReal & getSmoothingRadiusValue(){ return m_h; }
    const HReal & getSmoothingRadius() const { return m_h; }
    HReal & getTimeStepValue(){ return m_dt; }
    HReal & getMassValue(){ return m_mass; }
    HReal & getMeanDensityValue(){ return m_meanDensity;}
    HReal & getDensityFluctuationValue(){ return m_densityFluctuation;}
    HReal & getRealVolumeValue(){ return m_realVolume; }
    int & getParticleNumber(){ return m_particleNumber; }

    HReal & getViscosity(){return m_alpha;}
    const HReal & getViscosity() const {return m_alpha;}
    HReal & getFluidCohesion(){return m_fcohesion;}
    const HReal & getFluidCohesion() const {return m_fcohesion;}
    HReal & getBoundaryAdhesion(){return m_badhesion;}
    const HReal & getBoundaryAdhesion() const {return m_badhesion;}
    const HReal & getBoundaryFriction() const {return m_sigma;}
    HReal & getBoundaryFriction() {return m_sigma;}
    HReal & getTimeStep() {return m_dt;}
    const HReal& getTimeStep() const {return m_dt;}
};

bool pairCompare( const std::pair<int,int>& e1, const std::pair<int,int>& e2 );
}

#endif // SYSTEM_H
