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
#include "boundary.hpp"

#include "fluidParams.hpp"
#include "boundaryParams.hpp"

#include "gridUtility.hpp"
#include "triMesh.hpp"
#include "sampler.hpp"
#include "particleSource.hpp"
#include "kernel.hpp"

namespace hokusai
{

class System
{

public:
    System();
    System(int wishedParticleNumber);
    ~System();

public :

    //Stats
    int m_countTime;
    int m_countExport;
    HReal m_time;
    HReal m_meanDensity;
    HReal m_densityFluctuation;
    HReal m_realVolume;

    //Fluid parameters
    FluidParams m_fluidParams;
    BoundaryParams m_boundaryParams;

    //Simulation parameters
    int m_particleNumber;
    int m_boundaryNumber;

    Vec3r m_gravity;

    int m_maxPressureSolveIterationNb;
    HReal m_averageDensity;
    HReal m_dt;
    HReal m_maxDensityError;

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
    void computeWCSPHPressure(int i);
    void computeIISPHPressure(int i);

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


    //Initialize a dam break scenario
    const Vec3r& getGravity();
    void setGravity(const Vec3r& _gravity);
    void setParameters(int _number, HReal _volume=1.0, HReal _density=1000.0);

    void translateParticles(const Vec3r& t);
    void translateBoundaries(const Vec3r& t);

    //Source
    void addParticleSource(const ParticleSource& s);

    //Fluid sampling
    void addFluidParticle(const Vec3r& x, const Vec3r& v = Vec3r(0,0,0), const FluidParams& fluidParams=FluidParams());
    void addParticleBox(const Vec3r& offset, const Vec3r& dimension, const FluidParams& fluidParams=FluidParams());
    void addParticleSphere(const Vec3r& centre, const HReal radius, const FluidParams& fluidParams=FluidParams());

    //Boundary sampling
    void addBoundaryParticle(const Vec3r& x, const Vec3r& v = Vec3r(0,0,0), const BoundaryParams& boundaryParams=BoundaryParams());
    void addBoundaryMesh(const char* filename, const BoundaryParams& boundaryParams=BoundaryParams());
    void addBoundaryBox(const Vec3r& offset, const Vec3r& scale, const BoundaryParams& boundaryParams=BoundaryParams());
    void addBoundarySphere(const Vec3r& offset, const HReal& radius,const BoundaryParams& boundaryParams=BoundaryParams());
    void addBoundaryHemiSphere(const Vec3r& offset, const HReal& radius,const BoundaryParams& boundaryParams=BoundaryParams());
    void addBoundaryDisk(const Vec3r& offset, const HReal& radius,const BoundaryParams& boundaryParams=BoundaryParams());

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

    //Getter
    std::vector< Vec3r > getPosition();
    std::vector< Vec3r > getVelocity();
    std::vector< Vec3r > getNormal();
    std::vector< HReal > getDensity();
    std::vector< HReal > getMass();

    void write(const char* filename, std::vector< Vec3r > data);
    void write(const char* filename, std::vector<HReal> data);
    void exportState(const char* baseName);
    HReal getTime();
    HReal & getTimeStepValue();
    HReal & getMeanDensityValue();
    HReal & getDensityFluctuationValue();
    HReal & getRealVolumeValue();
    int & getParticleNumber();

    HReal & getTimeStep();
    const HReal& getTimeStep() const;
};

bool pairCompare( const std::pair<int,int>& e1, const std::pair<int,int>& e2 );
}

#endif // SYSTEM_H
