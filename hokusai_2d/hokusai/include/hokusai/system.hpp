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

#include <aljabr/Vec.hpp>
#include "utils.hpp"
#include "particle.hpp"
#include "gridUtility.hpp"
#include "triMesh.hpp"
#include "sampler.hpp"
#include "particleSource.hpp"

namespace hokusai
{

typedef aljabr::Vec3<double> Vec3r;
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

    double particlePerCell;
    double volume;
    double restDensity;
    double mean_density;
    double density_fluctuation;
    double real_volume;
    double mass;
    double h; // Smoothing radius
    double fcohesion; //Fluid cohesion
    double badhesion; //Boundary adhesion
    double sigma; //Boundary friction
    double cs;// Sound speed
    double alpha; // Viscosity
    double boundaryH;
    double dt;
    double time;
    double rho_avg_l;
    double maxEta;

    Vec gravity;

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
    void getNearestNeighbor(vector< int >& neighbors, const vector<vector<int> > &grid, const Vec& x);
    void getNearestNeighbor(const int i, const float radius);

    //Simulation Loop
    void prepareGrid();
    void computeSurfaceParticle();
    void predictAdvection();
    void predictRho(int i);
    void initializePressure(int i);
    void computeNormal(int i);
    bool isSurfaceParticle(int i, double treshold);
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
        Vec r = pi.x - pj.x;
        Vec vij = pi.v - pj.v;
        double dotVijRij = Vec::dotProduct(vij,r);
        if(dotVijRij < 0)
        {
            double kij = 2.0*restDensity/(pi.rho+pj.rho);
            double epsilon=0.01;
            Vec gradient(0.0);
            p_kernel.monaghanGradient(r, gradient);
            double Pij = -kij*(2.0*alpha*h*cs/(pi.rho+pj.rho)) * ( dotVijRij / (r.lengthSquared() + epsilon*h*h) );
            pi.f_adv += -kij*mass*mass*Pij*gradient;
        }
    }

    void computeBoundaryFrictionForces(int i, int j)
    {
        Particle& pi=particles[i];
            Boundary& bj=boundaries[j];
            Vec vij = pi.v;//-pj.v;
            Vec xij= pi.x - bj.x;
            double dotVijRij = Vec::dotProduct(vij,xij);
            if(dotVijRij<0)
            {
                Vec gradient(0.0);
                double epsilon=0.01;
                double nu = (sigma*h*cs)/(2.0*pi.rho);
                double Pij = -nu * ( std::min(dotVijRij,0.0) / (xij.lengthSquared() + epsilon*h*h) );
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
                Vec r = pi.x - pj.x;
                double kij = 2.0*restDensity/(pi.rho+pj.rho);
                double l = r.length();
                Vec cohesionForce = -(fcohesion*mass*mass*a_kernel.cohesionValue(l)/l) * r;
                Vec nij = pi.n-pj.n;
                Vec curvatureForce = -fcohesion*mass*nij;
                pi.f_adv += kij*(cohesionForce+curvatureForce);
            }
        }
    }
    
    void computeBoundaryAdhesionForces(int i, int j)
    {
        Particle& pi=particles[i];
            Boundary& bj=boundaries[j];
            Vec xij= pi.x - bj.x;
            double l = xij.length();
            pi.f_adv += -(badhesion*mass*boundaries[j].psi*a_kernel.adhesionValue(l)/l)*xij;
    }

    Vec computeDij(int i, int j)
    {
        Particle& pi=particles[i];
        Particle& pj=particles[j];
        Vec gradient(0.0);
        p_kernel.monaghanGradient(pi.x-pj.x, gradient);
        Vec d=-(dt*dt*mass)/pow(pj.rho,2)*gradient;
        return d;
    }

    void computePressure(int i);
    void computePressureForce(int i);

    void computeFluidPressureForce(int i, int j)
    {
        Vec gradient(0.0);
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
        Vec gradient(0.0);
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

    void addBoundaryParticle(const Vec& x, const Vec& v = Vec(0,0,0));
    void addFluidParticle(const Vec& x, const Vec& v = Vec(0,0,0));

    //Initialize a dam break scenario
    const Vec& getGravity();
    void setGravity(const Vec& _gravity);
    void setParameters(int _number, double _volume=1.0);
    void createParticleVolume(Vec& pos, double width, double height, double depth, double spacing, int particleMax);

    void translateParticles(const Vec& t);
    void translateBoundaries(const Vec& t);

    void addParticleBox(double width, double height, double depth, double spacing);
    void addParticleBox(const Vec& offset, const Vec& dimension);
    void addParticleSphere(const Vec& centre, const double radius);

    void addParticleSource(const ParticleSource& s);

    //Boundary sampling
    void addBoundaryMesh(const char* filename);
    void addBoundaryBox(const Vec& offset, const Vec& scale);
    void addBoundarySphere(const Vec& offset, const double& radius);
    void addBoundaryHemiSphere(const Vec& offset, const double& radius);
    void addBoundaryDisk(const Vec& offset, const double& radius);

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
    vector< double > getDensity(){ vector<double> density; for(int i=0; i<particleNumber; ++i){density.push_back(particles[i].rho);} return density;}
    vector< double > getMass(){ vector<double> o_mass; for(int i=0; i<particleNumber; ++i){o_mass.push_back(mass);} return o_mass;}

    void write(const char* filename, vector< Vec3r > data);
    void write(const char* filename, vector<double> data);
    void exportState(const char* baseName);
    double getTime(){return time;}
    double & getSmoothingRadiusValue(){ return h; }
    const double & getSmoothingRadius() const { return h; }
    double & getTimeStepValue(){ return dt; }
    double & getMassValue(){ return mass; }
    double & getMeanDensityValue(){ return mean_density;}
    double & getDensityFluctuationValue(){ return density_fluctuation;}
    double & getRealVolumeValue(){ return real_volume; }
    int & getParticleNumber(){ return particleNumber; }

    SReal & getViscosity(){return alpha;}
    const SReal & getViscosity() const {return alpha;}
    SReal & getFluidCohesion(){return fcohesion;}
    const SReal & getFluidCohesion() const {return fcohesion;}
    SReal & getBoundaryAdhesion(){return badhesion;}
    const SReal & getBoundaryAdhesion() const {return badhesion;}
    const SReal & getBoundaryFriction() const {return sigma;}
    SReal & getBoundaryFriction() {return sigma;}
    SReal & getTimeStep() {return dt;}
    const SReal& getTimeStep() const {return dt;}

    void applyShepardFilter();
};

bool pairCompare( const std::pair<int,int>& e1, const std::pair<int,int>& e2 );
}

#endif // SYSTEM_H