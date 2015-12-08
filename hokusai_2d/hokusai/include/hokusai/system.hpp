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

#include "utility.hpp"
#include "gridUtility.hpp"
#include "utils.hpp"
#include "particle.hpp"
#include "particleSource.hpp"
#include "particleSink.hpp"

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

    Vec2d gravity;

    AkinciKernel a_kernel;
    MonaghanKernel p_kernel;
    BoundaryKernel b_kernel;

    vector<Particle> particles;
    vector<Boundary> boundaries;
    vector<MovingBoundary> m_moving_boundaries;

    Grid2dUtility gridInfo;
    vector< vector<int> > boundaryGrid;
    vector< vector<int> > fluidGrid;
    std::vector<ParticleSink> p_sinks;
    std::vector<ParticleSource> p_sources;

public :
    //Simulation Loop

    ///Compute fluid particles which are in a given radius from a given position.
    void getNearestFluidNeighbor(vector< int >& neighbors, const Vec2d& x, const double radius);

    ///Compute fluid and boundary neighbors for a particle i in a given radius.
    void getNearestNeighbor(const int i, const double radius);

    ///Compute neighbor for a position x in a given radius from a given grid.
    void getNearestNeighbor(vector< int >& neighbor, const vector< vector<int> >& grid, const Vec2d &x);

    ///Clean fluid/boundary grid, resize and fill them
    void prepareGrid();

    void predictAdvection();
    void predictRho(int i);
    void initializePressure(int i);
    void computeNormal(int i);
    void computeRho(int i);
    void computeAdvectionForces(int i);
    void predictVelocity(int i);
    void computeDii(int i);
    void computeAii(int i);
    void pressureSolve();
    void computeSumDijPj(int i);

    void addParticleSink(const ParticleSink& s);
    void addParticleSource(const ParticleSource& s);
    void neighborNumberSurfaceDetection(std::vector<int>& surface, int threshold);
    void massCenterDistanceSurfaceDetection(std::vector<int>& surface, double distanceThreshold, int neighborNumberThreshold);
    void colorFieldSurfaceDetection(std::vector<int>& surface, double normalThreshold, int neighborNumberThreshold);

    void set_benchmark_sink(double timeStep);
    void set_benchmark_source(double timeStep);
    void set_benchmark_dambreak(double timeStep);
    void set_benchmark_wave(double timeStep);
    void set_benchmark_splash(double timeStep);

    inline Vec2d computeVelocity(const Vec2d& position, const std::vector<int>& neighbors)
    {
        Vec2d velocity(0.0,0.0);
        for(const int& i : neighbors)
        {
            Particle& pi=particles[i];
            velocity += pi.v*(mass/pi.rho)*p_kernel.monaghanValue(position-pi.x);
        }
        return velocity;
    }

    inline void computeViscosityForces(int i, int j)
    {
        Particle& pi=particles[i];
        Particle& pj=particles[j];
        Vec2d r = pi.x - pj.x;
        Vec2d vij = pi.v - pj.v;
        double dotVijRij = vij.dot(r);
        if(dotVijRij < 0)
        {
            double kij = 2.0*restDensity/(pi.rho+pj.rho);
            double epsilon=0.01;
            Vec2d gradient(0.0,0.0);
            p_kernel.monaghanGradient(r, gradient);
            double Pij = -kij*(2.0*alpha*h*cs/(pi.rho+pj.rho)) * ( dotVijRij / (r.squaredNorm() + epsilon*h*h) );
            pi.f_adv += -kij*mass*mass*Pij*gradient;
        }
    }

    inline void computeSurfaceTensionForces(int i, int j)
    {
        if(i!=j)
        {
            Particle& pi=particles[i];
            Particle& pj=particles[j];
            Vec2d r = pi.x - pj.x;
            double kij = 2.0*restDensity/(pi.rho+pj.rho);
            double l = r.norm();
            Vec2d cohesionForce = -(fcohesion*mass*mass*a_kernel.cohesionValue(l)/l) * r;
            Vec2d nij = pi.n-pj.n;
            Vec2d curvatureForce = -fcohesion*mass*nij;
            pi.f_adv += kij*(cohesionForce+curvatureForce);
        }
    }

    inline void computeBoundaryFrictionForces(int i, int j)
    {
        Particle& pi=particles[i];
        Boundary& bj=boundaries[j];
        Vec2d vij = pi.v;//-pj.v;
        Vec2d xij= pi.x - bj.x;
        double dotVijRij = vij.dot(xij);
        if(dotVijRij<0)
        {
            Vec2d gradient(0.0,0.0);
            double epsilon=0.01;
            double nu = (sigma*h*cs)/(2.0*pi.rho);
            double Pij = -nu * ( std::min(dotVijRij,0.0) / (xij.squaredNorm() + epsilon*h*h) );
            p_kernel.monaghanGradient(xij, gradient);
            pi.f_adv += -mass*bj.psi*Pij*gradient;
        }
    }
    
    inline void computeBoundaryAdhesionForces(int i, int j)
    {
        Particle& pi=particles[i];
        Boundary& bj=boundaries[j];
        Vec2d xij= pi.x - bj.x;
        double l = xij.norm();
        pi.f_adv += -(badhesion*mass*boundaries[j].psi*a_kernel.adhesionValue(l)/l)*xij;
    }

    inline Vec2d computeDij(int i, int j)
    {
        Particle& pi=particles[i];
        Particle& pj=particles[j];
        Vec2d gradient(0.0,0.0);
        p_kernel.monaghanGradient(pi.x-pj.x, gradient);
        Vec2d d=-(dt*dt*mass)/pow(pj.rho,2)*gradient;
        return d;
    }

    void computePressure(int i);
    void computePressureForce(int i);

    inline void computeFluidPressureForce(int i, int j)
    {
        Vec2d gradient(0.0,0.0);
        Particle& pi=particles[i];
        Particle& pj=particles[j];
        p_kernel.monaghanGradient(pi.x-pj.x, gradient);
        if( i!=j )
        {
            pi.f_p += -mass*mass*( pi.p/pow(pi.rho,2) + pj.p/pow(pj.rho,2) ) * gradient;
        }
    }

    inline void computeBoundaryPressureForce(int i, int j)
    {
        Vec2d gradient(0.0,0.0);
        Particle& pi=particles[i];
        Boundary& bj=boundaries[j];
        p_kernel.monaghanGradient(pi.x-bj.x, gradient);
        pi.f_p += -mass*bj.psi*( pi.p/pow(pi.rho,2) ) * gradient;
    }

    void updateGridInfo();
    void updateGridInfo(Vec2d& minBB, Vec2d& maxBB);
    void computeBoundingBox(Vec2d &minBB, Vec2d &maxBB);

    void computeError();
    void integration();
    void simulate();
    void applySources();
    void applySinks();
    void applyMovingBoundaries();

    void addBoundaryParticle(const Vec2d& x, const Vec2d& v);
    void addFluidParticle(const Vec2d& x, const Vec2d& v);

    //Initialize a dam break scenario
    Vec2d& getGravity(){return gravity;}
    const Vec2d& getGravity() const;
    void setGravity(const Vec2d& _gravity);
    void setParameters(int _number, double _volume=1.0);
    void updateParameters();
    void createParticleVolume(Vec2d& pos, double width, double height, double depth, double spacing, int particleMax);

    void translateParticles(const Vec2d& t);
    void translateBoundaries(const Vec2d& t);

    void cleanFluidParticle();

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
    vector< Vec2d > getBoundaryPosition(){ vector<Vec2d > pos; for(int i=0; i<boundaryNumber; ++i){pos.push_back(boundaries[i].x);} return pos;}
    vector< Vec2d > getFluidPosition(){ vector<Vec2d > pos; for(int i=0; i<particleNumber; ++i){pos.push_back(particles[i].x);} return pos;}
    vector< Vec2d > getVelocity(){ vector<Vec2d > vel; for(int i=0; i<particleNumber; ++i){vel.push_back(particles[i].v);} return vel;}
    vector< Vec2d > getNormal(){ vector<Vec2d > normal; for(int i=0; i<particleNumber; ++i){normal.push_back(particles[i].n);} return normal;}
    vector< double > getDensity(){ vector<double> density; for(int i=0; i<particleNumber; ++i){density.push_back(particles[i].rho);} return density;}
    vector< double > getMass(){ vector<double> o_mass; for(int i=0; i<particleNumber; ++i){o_mass.push_back(mass);} return o_mass;}

    void write(const char* filename, vector< Vec2d > data);
    void write(const char* filename, vector<double> data);
    void exportState(const char* baseName);
    double getTime(){return time;}

    double & getSmoothingRadius(){ return h; }
    const double & getSmoothingRadius() const { return h; }
    void setSmoothingRadius( double _h){ h = _h; }

    double & getTimeStep() {return dt;}
    const double& getTimeStep() const {return dt;}
    void setTimeStep(double _dt){dt = _dt;}

    double & getRestDensity(){ return restDensity;}
    const double & getRestDensity() const{ return restDensity;}
    void setRestDensity( const double _restDensity) { restDensity = _restDensity;}

    double & getMassValue(){ return mass; }
    double & getMeanDensityValue(){ return mean_density;}
    double & getDensityFluctuationValue(){ return density_fluctuation;}
    double & getRealVolumeValue(){ return real_volume; }
    int & getParticleNumber(){ return particleNumber; }

    double & getViscosity(){return alpha;}
    const double & getViscosity() const {return alpha;}
    void setViscosity(double _alpha){alpha = _alpha;}

    double & getFluidCohesion(){return fcohesion;}
    const double & getFluidCohesion() const {return fcohesion;}
    void setFluidCohesion(double _fcohesion){fcohesion = _fcohesion;}

    double & getBoundaryAdhesion(){return badhesion;}
    const double & getBoundaryAdhesion() const {return badhesion;}
    void setBoundaryAdhesion(double _badhesion){badhesion = _badhesion;}

    const double & getBoundaryFriction() const {return sigma;}
    double & getBoundaryFriction() {return sigma;}
    void setBoundaryFriction(double _sigma){sigma = _sigma;}

    void addMovingBoundaryBox(Vec2d& offset, Vec2d& scale, Vec2d &translation, bool& oscillatory, Vec2d& delta);
    std::vector<int> addBoundaryBox(Vec2d &offset, Vec2d &scale);
    void addParticleBox(Vec2d offset, Vec2d scale);
};

}

#endif // SYSTEM_H
