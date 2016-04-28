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
#include "iisph.hpp"
#include "kernel.hpp"

namespace hokusai
{

class System
{

public:

    System();
    System(int resolution);
    ~System();

public :

    //Stats
    int m_countTime;
    int m_countExport;
    HReal m_time;
    HReal m_volume;
    HReal m_restDensity;
    HReal m_meanDensity;
    HReal m_densityFluctuation;
    HReal m_realVolume;

    //Fluid parameters
    HReal m_mass;
    HReal m_h; // Smoothing radius
    HReal m_fcohesion; //Fluid cohesion
    HReal m_badhesion; //Boundary adhesion
    HReal m_sigma; //Boundary friction
    HReal m_cs;// Sound speed
    HReal m_alpha; // Viscosity
    HReal m_boundaryH;

    //Simulation parameters
    int m_particleNumber;
    int m_boundaryNumber;
    HReal m_particlePerCell;
    HReal m_dt;

    Vec3r m_gravity;

    IISPHParams m_iisphParams;

    std::vector<Particle> m_particles;
    std::vector<Boundary> m_boundaries;

    GridUtility m_gridInfo;
    std::vector< std::vector<int> > m_boundaryGrid;
    std::vector< std::vector<int> > m_fluidGrid;

    std::vector< ParticleSource > m_pSources;

public :
    void getNearestNeighbor(std::vector< int >& neighbors, const std::vector<std::vector<int> > &grid, const Vec3r& x);
    void getNearestNeighbor(const int i, const HReal radius);

    //Simulation Loop
    void prepareGrid();
    void computeSurfaceParticle();
    void predictAdvection();
    void predictDensity(int i);
    void initializePressure(int i);
    void computeNormal(int i);
    bool isSurfaceParticle(int i, HReal treshold);
    std::vector<Particle> getSurfaceParticle();
    void computeDensity(int i);
    void predictVelocity(int i);
    void computeDii(int i);
    void computeDii_Fluid(int i);
    void computeDii_Boundary(int i);
    void computeAii(int i);
    void pressureSolve();
    void computeSumDijPj(int i);    
    Vec3r computeDij(int i, int j);
    void computePressure(int i);

    void computeAdvectionForces(int i);
    void computeViscosityForces(int i, int j);
    void computeBoundaryFrictionForces(int i, int j);
    void computeSurfaceTensionForces(int i, int j);
    void computeBoundaryAdhesionForces(int i, int j);
    void computePressureForce(int i);
    void computeFluidPressureForce(int i, int j);
    void computeBoundaryPressureForce(int i, int j);

    void computeError();
    void integration();
    void computeSimulationStep();
    void applySources();
    void applySinks();

    void addBoundaryParticle(const Vec3r& x, const Vec3r& v = Vec3r(0,0,0));
    void addFluidParticle(const Vec3r& x, const Vec3r& v = Vec3r(0,0,0));

    //Initialize a dam break scenario
    const Vec3r& getGravity();
    void setGravity(const Vec3r& _gravity);
    void setParameters(int _number, HReal _volume=1.0);

    void translateParticles(const Vec3r& t);
    void translateBoundaries(const Vec3r& t);

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
    std::vector< Vec3r > getPosition(){ std::vector<Vec3r > pos; for(int i=0; i<m_particleNumber; ++i){pos.push_back(m_particles[i].x);} return pos;}
    std::vector< Vec3r > getVelocity(){ std::vector<Vec3r > vel; for(int i=0; i<m_particleNumber; ++i){vel.push_back(m_particles[i].v);} return vel;}
    std::vector< Vec3r > getNormal(){ std::vector<Vec3r > normal; for(int i=0; i<m_particleNumber; ++i){normal.push_back(m_particles[i].n);} return normal;}
    std::vector< HReal > getDensity(){ std::vector<HReal> density; for(int i=0; i<m_particleNumber; ++i){density.push_back(m_particles[i].rho);} return density;}
    std::vector< HReal > getMass(){ std::vector<HReal> o_mass; for(int i=0; i<m_particleNumber; ++i){o_mass.push_back(m_mass);} return o_mass;}

    void write(const char* filename, std::vector< Vec3r > data);
    void write(const char* filename, std::vector<HReal> data);
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
