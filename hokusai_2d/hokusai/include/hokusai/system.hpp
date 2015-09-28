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
#include "utils.hpp"
#include "particle.hpp"

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

    //AkinciKernel a_kernel;
    MonaghanKernel p_kernel;
    BoundaryKernel b_kernel;

    vector<Particle> particles;
    vector<Boundary> boundaries;

public :
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
        Vec2d r = pi.x - pj.x;
        Vec2d vij = pi.v - pj.v;
        double dotVijRij = vij.dot(r);
        if(dotVijRij < 0)
        {
            double kij = 2.0*restDensity/(pi.rho+pj.rho);
            double epsilon=0.01;
            Vec2d gradient(0.0);
            p_kernel.monaghanGradient(r, gradient);
            double Pij = -kij*(2.0*alpha*h*cs/(pi.rho+pj.rho)) * ( dotVijRij / (r.squaredNorm() + epsilon*h*h) );
            pi.f_adv += -kij*mass*mass*Pij*gradient;
        }
    }

    void computeBoundaryFrictionForces(int i, int j)
    {
        Particle& pi=particles[i];
            Boundary& bj=boundaries[j];
            Vec2d vij = pi.v;//-pj.v;
            Vec2d xij= pi.x - bj.x;
            double dotVijRij = vij.dot(xij);
            if(dotVijRij<0)
            {
                Vec2d gradient(0.0);
                double epsilon=0.01;
                double nu = (sigma*h*cs)/(2.0*pi.rho);
                double Pij = -nu * ( std::min(dotVijRij,0.0) / (xij.squaredNorm() + epsilon*h*h) );
                p_kernel.monaghanGradient(xij, gradient);
                pi.f_adv += -mass*bj.psi*Pij*gradient;
            }
    }

    void computeSurfaceTensionForces(int i, int j)
    {
//        if(i!=j)
//        {
//            Particle& pi=particles[i];
//            Particle& pj=particles[j];
//            if(pi.isSurface==true || pj.isSurface==true)
//            {
//                Vec2d r = pi.x - pj.x;
//                double kij = 2.0*restDensity/(pi.rho+pj.rho);
//                double l = r.norm();
//                Vec2d cohesionForce = -(fcohesion*mass*mass*a_kernel.cohesionValue(l)/l) * r;
//                Vec2d nij = pi.n-pj.n;
//                Vec2d curvatureForce = -fcohesion*mass*nij;
//                pi.f_adv += kij*(cohesionForce+curvatureForce);
//            }
//        }
    }
    
    void computeBoundaryAdhesionForces(int i, int j)
    {
//        Particle& pi=particles[i];
//            Boundary& bj=boundaries[j];
//            Vec2d xij= pi.x - bj.x;
//            double l = xij.norm();
//            pi.f_adv += -(badhesion*mass*boundaries[j].psi*a_kernel.adhesionValue(l)/l)*xij;
    }

    Vec2d computeDij(int i, int j)
    {
        Particle& pi=particles[i];
        Particle& pj=particles[j];
        Vec2d gradient(0.0);
        p_kernel.monaghanGradient(pi.x-pj.x, gradient);
        Vec2d d=-(dt*dt*mass)/pow(pj.rho,2)*gradient;
        return d;
    }

    void computePressure(int i);
    void computePressureForce(int i);

    void computeFluidPressureForce(int i, int j)
    {
        Vec2d gradient(0.0);
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
        Vec2d gradient(0.0);
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

    void addBoundaryParticle(const Vec2d& x, const Vec2d& v);
    void addFluidParticle(const Vec2d& x, const Vec2d& v);

    //Initialize a dam break scenario
    const Vec2d& getGravity();
    void setGravity(const Vec2d& _gravity);
    void setParameters(int _number, double _volume=1.0);
    void createParticleVolume(Vec2d& pos, double width, double height, double depth, double spacing, int particleMax);

    void translateParticles(const Vec2d& t);
    void translateBoundaries(const Vec2d& t);

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
    double & getSmoothingRadiusValue(){ return h; }
    const double & getSmoothingRadius() const { return h; }
    double & getTimeStepValue(){ return dt; }
    double & getMassValue(){ return mass; }
    double & getMeanDensityValue(){ return mean_density;}
    double & getDensityFluctuationValue(){ return density_fluctuation;}
    double & getRealVolumeValue(){ return real_volume; }
    int & getParticleNumber(){ return particleNumber; }

    double & getViscosity(){return alpha;}
    const double & getViscosity() const {return alpha;}
    double & getFluidCohesion(){return fcohesion;}
    const double & getFluidCohesion() const {return fcohesion;}
    double & getBoundaryAdhesion(){return badhesion;}
    const double & getBoundaryAdhesion() const {return badhesion;}
    const double & getBoundaryFriction() const {return sigma;}
    double & getBoundaryFriction() {return sigma;}
    double & getTimeStep() {return dt;}
    const double& getTimeStep() const {return dt;}

    void addBoundaryBox(Vec2d offset, Vec2d scale);
    void addParticleBox(Vec2d offset, Vec2d scale);
};

}

#endif // SYSTEM_H
