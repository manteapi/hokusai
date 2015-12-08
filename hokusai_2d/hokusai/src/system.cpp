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

#include "../include/hokusai/system.hpp"
#include "../include/hokusai/utils.hpp"
#include <algorithm>
#include <assert.h>

#ifdef HOKUSAI_USING_OPENMP
#include <omp.h>
#endif

#include <fstream>
#include <sstream>

namespace hokusai
{
//-----------------------------------------------------
//         System : creation and simulation functions
//-----------------------------------------------------

System::System()
{
    p_sinks.clear();
    p_sources.clear();
    countExport = 0;
    countTime = 0;
    particleNumber = 0;
    boundaryNumber = 0;

    volume = 0.0;
    restDensity = 0.0;
    mean_density = 0.0;
    density_fluctuation = 0.0;
    real_volume = 0.0;
    mass = 0.0;
    h = 0.0;
    fcohesion = 0.0;
    badhesion = 0.0;
    cs = 0.0;
    alpha = 0.0;
    boundaryH = 0.0;
    dt = 0.0;
    time = 0.0;
    rho_avg_l = 0.0;
    maxEta = 1.0;

    gravity = Vec2d(0,-9.81);

    a_kernel = AkinciKernel();
    p_kernel = MonaghanKernel();
    b_kernel = BoundaryKernel();

    gridInfo = Grid2dUtility();
    particles = vector<Particle>();
    boundaries = vector<Boundary>();
    m_moving_boundaries = vector<MovingBoundary>();
}

System::System(int resolution)
{
    p_sinks.clear();
    p_sources.clear();
    countExport = 0;
    countTime = 0;
    particleNumber = 0;
    boundaryNumber = 0;

    volume = 0.0;
    restDensity = 0.0;
    mean_density = 0.0;
    density_fluctuation = 0.0;
    real_volume = 0.0;
    mass = 0.0;
    h = 0.0;
    fcohesion = 0.0;
    badhesion = 0.0;
    cs = 0.0;
    alpha = 0.0;
    boundaryH = 0.0;
    dt = 0.0;
    time = 0.0;
    rho_avg_l = 0.0;
    maxEta = 1.0;

    gravity = Vec2d(0,-9.81);

    a_kernel = AkinciKernel();
    p_kernel = MonaghanKernel();
    b_kernel = BoundaryKernel();

    gridInfo = Grid2dUtility();
    particles = vector<Particle>();
    boundaries = vector<Boundary>();
    m_moving_boundaries = vector<MovingBoundary>();

    setParameters(resolution, 1.0);
}

System::~System(){}

void System::neighborNumberSurfaceDetection(std::vector<int>& surface, int threshold)
{
    for(size_t i=0; i<particles.size(); ++i)
    {
        Particle& pi = particles[i];
        int neighborNumber=0;
        for(size_t j=0; j<pi.fluidNeighbor.size(); ++j)
        {
            Particle& pj = particles[pi.fluidNeighbor[j]];
            if(p_kernel.monaghanValue(pi.x-pj.x)>0)
            {
                neighborNumber++;
            }
        }
        if(neighborNumber<=threshold)
        {
            surface.push_back(i);
        }
    }
}

void System::massCenterDistanceSurfaceDetection(std::vector<int>& surface, double distanceThreshold, int neighborNumberThreshold)
{
    for(size_t i=0; i<particles.size(); ++i)
    {
        Particle& pi = particles[i];
	    Vec2d x = pi.x;
	    Vec2d massCenter(0.0,0.0);
	    int neighborNumber = 0;
        for(size_t j=0; j<pi.fluidNeighbor.size(); ++j)
        {
	        int neighborId = pi.fluidNeighbor[j];
            Particle& pj = particles[neighborId];
            if(p_kernel.monaghanValue(pi.x-pj.x)>0)
            {
                neighborNumber++;
	            massCenter += particles[neighborId].x;
            }
        }
	    massCenter/=neighborNumber;
	    double distance = (x-massCenter).norm();
	    if(distance>=distanceThreshold || neighborNumber<=neighborNumberThreshold)
            surface.push_back(i);
    }
}
    
void System::colorFieldSurfaceDetection(std::vector<int>& surface, double normalThreshold, int neighborNumberThreshold)
{
    for(size_t i=0; i<particles.size(); ++i)
    {
	    Particle& pi = particles[i];
	    int neighborNumber = 0;
        double normalLength = pi.n.norm();
        for(size_t j=0; j<pi.fluidNeighbor.size(); ++j)
        {
	        int neighborId = pi.fluidNeighbor[j];
            Particle& pj = particles[neighborId];
            if(p_kernel.monaghanValue(pi.x-pj.x)>0)
            {
                neighborNumber++;
            }
        }
        if(neighborNumber<=neighborNumberThreshold || normalLength>=normalThreshold)
        {
            surface.push_back(i);
        }
    }
}

void System::addParticleSink(const ParticleSink& s)
{
    p_sinks.push_back(s);
}

void System::addParticleSource(const ParticleSource& s)
{
    p_sources.push_back(s);
}

void System::applySources()
{
    for(ParticleSource& s : p_sources)
    {
        std::vector<Particle> p_new = s.apply(this->time);
        for(const Particle& p : p_new)
        {
            particles.push_back(p);
            particleNumber++;
        }
    }
}

void System::applySinks()
{
    std::vector<int> p_remove;
    for(ParticleSink& s : p_sinks)
    {
        Vec2i minBB = gridInfo.worldToGrid(s.m_minBB);
        Vec2i maxBB = gridInfo.worldToGrid(s.m_maxBB);
        for(int i=minBB[0]; i<=maxBB[0]; ++i)
        {
            for(int j=minBB[1]; j<=maxBB[1]; ++j)
            {
                int cellId = gridInfo.cellId(i,j);
                if(gridInfo.isInside(cellId))
                {
                    for(size_t k=0; k<fluidGrid[cellId].size(); ++k)
                    {
                        int pId = fluidGrid[cellId][k];
                        if(s.contain(particles[pId].x))
                        {
                            p_remove.push_back(pId);
                        }
                    }
                }
            }
        }
    }
    for(const int& pId : p_remove)
    {
        particles.erase(particles.begin()+pId);
        particleNumber--;
    }
}

void System::computeRho(int i)
{
    Particle& pi=particles[i];
    vector<int>& fneighbors=pi.fluidNeighbor;
    vector<int>& bneighbors=pi.boundaryNeighbor;
    pi.rho=0.0;
    for(int& j : fneighbors)
    {
        Particle& pj=particles[j];
        pi.rho += mass*p_kernel.monaghanValue(pi.x-pj.x);
    }
    for(int& j : bneighbors)
    {
        Boundary& bj = boundaries[j];
        pi.rho += p_kernel.monaghanValue(pi.x-bj.x)*bj.psi;
    }
}

void System::computeNormal(int i)
{
    //Compute normal
    Particle& pi=particles[i];
    vector< int > & neighbors = particles[i].fluidNeighbor;
    Vec2d n(0.0,0.0);
    Vec2d gradient(0.0,0.0);
    for(int& j : neighbors)
    {
        if(i!=j)
        {
            Particle& pj = particles[j];
            p_kernel.monaghanGradient(pi.x-pj.x, gradient);
            n += (mass/pj.rho)*gradient;
        }
    }
    pi.n = h*n;
}

void System::computeAdvectionForces(int i)
{
    Particle& pi=particles[i];
    pi.f_adv = Vec2d::Zero();

    Vec2d r(0.0,0.0), vij(0.0,0.0), gradient(0.0,0.0), cohesionForce(0.0,0.0), curvatureForce(0.0,0.0), nij(0.0,0.0);
    double dotVijRij=0.0, kij=0.0, l=0.0, Pij=0.0, epsilon=0.01, nu=0.0;

    ///Fluid forces : Viscosity + surface tension
    for(int& j : pi.fluidNeighbor)
    {
        Particle& pi=particles[i];
        Particle& pj=particles[j];
        r = pi.x - pj.x;
        vij = pi.v - pj.v;
        dotVijRij = vij.dot(r);
        kij = 2.0*restDensity/(pi.rho+pj.rho);
        l = r.norm();

        ///Viscosity forces
        if(dotVijRij < 0)
        {
            p_kernel.monaghanGradient(r, gradient);
            Pij = -kij*(2.0*alpha*h*cs/(pi.rho+pj.rho)) * ( dotVijRij / (l*l + epsilon*h*h) );
            pi.f_adv += -kij*mass*mass*Pij*gradient;
        }

        ///Surface tension forces
        if(i!=j)
        {
            cohesionForce = -(fcohesion*mass*mass*a_kernel.cohesionValue(l)/l) * r;
            nij = pi.n-pj.n;
            curvatureForce = -fcohesion*mass*nij;
            pi.f_adv += kij*(cohesionForce+curvatureForce);
        }
    }

    ///Boundary forces : Friction + adhesion
    for(int& j : pi.boundaryNeighbor)
    {
        Particle& pi=particles[i];
        Boundary& bj=boundaries[j];
        vij = pi.v;//-pj.v;
        r= pi.x - bj.x;
        l = r.norm();
        dotVijRij = vij.dot(r);

        ///Friction
        if(dotVijRij<0)
        {
            Vec2d gradient(0.0,0.0);
            nu = (sigma*h*cs)/(2.0*pi.rho);
            Pij = -nu * ( std::min(dotVijRij,0.0) / (l*l + epsilon*h*h) );
            p_kernel.monaghanGradient(r, gradient);
            pi.f_adv += -mass*bj.psi*Pij*gradient;
        }

        ///Adhesion
        pi.f_adv += -(badhesion*mass*boundaries[j].psi*a_kernel.adhesionValue(l)/l)*r;
    }

    pi.f_adv+=gravity*mass;
}

void System::predictVelocity(int i)
{
    Particle& pi=particles[i];
    pi.v_adv = pi.v + (dt/mass)*pi.f_adv;
}

void System::predictRho(int i)
{
    Particle& pi=particles[i];
    vector<int>& fneighbors=pi.fluidNeighbor;
    vector<int>& bneighbors=pi.boundaryNeighbor;
    double fdrho=0.0, bdrho=0.0;
    Vec2d gradient(0.0,0.0), vij_adv(0.0,0.0), vb(0.1,0.1), v(0.0,0.0);

    for(int& j : fneighbors)
    {
        if(i!=j)
        {
            Particle& pj=particles[j];
            p_kernel.monaghanGradient(pi.x-pj.x, gradient);
            vij_adv=pi.v_adv-pj.v_adv;
            fdrho+=mass*vij_adv.dot(gradient);
        }
    }

    for(int& j: bneighbors)
    {
        Boundary& bj=boundaries[j];
        v = pi.v_adv - vb; //vb(t+dt)
        p_kernel.monaghanGradient(pi.x-bj.x, gradient);
        bdrho+=bj.psi*v.dot(gradient);
    }

    pi.rho_adv = pi.rho + dt*( fdrho + bdrho );
}

void System::computeSumDijPj(int i)
{
    Particle& pi=particles[i];
    vector<int>& fneighbors=pi.fluidNeighbor;
    pi.sum_dij = Vec2d::Zero();
    Vec2d gradient(0.0,0.0);

    for(int& j : fneighbors)
    {
        if(i!=j)
        {
            Particle& pj=particles[j];
            p_kernel.monaghanGradient(pi.x-pj.x, gradient);
            pi.sum_dij+=(-mass/pow(pj.rho,2))*pj.p_l*gradient;
        }
    }
    pi.sum_dij *= pow(dt,2);
}

void System::computePressure(int i)
{
    Particle& pi=particles[i];
    vector<int>& fneighbors=pi.fluidNeighbor, bneighbors=pi.boundaryNeighbor;
    double fsum=0.0, bsum=0.0, omega=0.5;
    Vec2d gradient_ij(0.0,0.0), dji(0.0,0.0), aux(0.0,0.0);

    for(int& j : fneighbors)
    {
        if(i!=j)
        {
            Particle& pj=particles[j];
            dji=computeDij(j, i);
            p_kernel.monaghanGradient(pi.x-pj.x, gradient_ij);
            aux = pi.sum_dij - (pj.dii_fluid+pj.dii_boundary)*pj.p_l - (pj.sum_dij - dji*pi.p_l);
            fsum+=mass*aux.dot(gradient_ij);
        }
    }

    for(int& j : bneighbors)
    {
        Boundary& bj=boundaries[j];
        p_kernel.monaghanGradient(pi.x-bj.x, gradient_ij);
        bsum+=bj.psi*pi.sum_dij.dot(gradient_ij);
    }

    double previousPl = pi.p_l;
    pi.rho_corr = pi.rho_adv + fsum + bsum;
    if(std::abs(pi.aii)>std::numeric_limits<double>::epsilon())
        pi.p_l = (1-omega)*previousPl + (omega/pi.aii)*(restDensity - pi.rho_corr);
    else
        pi.p_l = 0.0;
    pi.p = std::max(pi.p_l,0.0);
    pi.p_l = pi.p;
    pi.rho_corr += pi.aii*previousPl;
}

void System::computePressureForce(int i)
{
    Particle& pi=particles[i];
    pi.f_p = Vec2d::Zero();
    Vec2d gradient(0.0,0.0);

    //Fluid Pressure Force
    for(int & j : particles[i].fluidNeighbor)
    {
        Particle& pi=particles[i];
        Particle& pj=particles[j];
        p_kernel.monaghanGradient(pi.x-pj.x, gradient);
        if( i!=j )
        {
            pi.f_p += -mass*mass*( pi.p/pow(pi.rho,2) + pj.p/pow(pj.rho,2) ) * gradient;
        }
    }

    //Boundary Pressure Force [Akinci 2012]
    for(int & j : particles[i].boundaryNeighbor)
    {
        Particle& pi=particles[i];
        Boundary& bj=boundaries[j];
        p_kernel.monaghanGradient(pi.x-bj.x, gradient);
        pi.f_p += -mass*bj.psi*( pi.p/pow(pi.rho,2) ) * gradient;
    }
}

void System::initializePressure(int i)
{
    Particle& pi=particles[i];
    pi.p_l=0.5*pi.p;
}

void System::computeError()
{
    rho_avg_l=0.0;
    for(int i=0; i<particleNumber; ++i)
        rho_avg_l+=particles[i].rho_corr;
    rho_avg_l /= particleNumber;
}

void System::updateGridInfo()
{
    Vec2d minBB, maxBB;
    computeBoundingBox(minBB, maxBB);
    updateGridInfo(minBB, maxBB);
}

void System::updateGridInfo (Vec2d &minBB, Vec2d &maxBB)
{
    double eOffset = 2.0*h;
    double eScale = 4.0*h;
    Vec2d gOffset = minBB-Vec2d(eOffset,eOffset);
    Vec2d gScale = (maxBB-minBB)+Vec2d(eScale,eScale);
    gridInfo.update(gOffset, gScale, 2.0*h);
}

void System::computeBoundingBox(Vec2d& minBB, Vec2d& maxBB)
{
    for(size_t i=0; i<2; ++i)
    {
        minBB[i] = std::numeric_limits<double>::max();
        maxBB[i] = -std::numeric_limits<double>::max();
    }
    for(Boundary& pi : boundaries)
    {
        for(size_t i=0; i<2; ++i)
        {
            minBB[i] = (pi.x[i]<minBB[i]) ? pi.x[i]  : minBB[i];
            maxBB[i] = (pi.x[i]>maxBB[i]) ? pi.x[i]  : maxBB[i];
        }
    }
}

void System::computeDii(int i)
{
    Particle& pi=particles[i];
    pi.dii_fluid = Vec2d::Zero();
    pi.dii_boundary = Vec2d::Zero();
    Vec2d gradient(0.0,0.0);

    for(int& j : pi.fluidNeighbor)
    {
        if(i!=j)
        {
            Particle& pj=particles[j];
            p_kernel.monaghanGradient(pi.x-pj.x, gradient);
            pi.dii_fluid+=(-dt*dt*mass/pow(pi.rho,2))*gradient;
        }
    }

    for(int& j : pi.boundaryNeighbor)
    {
        Boundary& bj=boundaries[j];
        p_kernel.monaghanGradient(pi.x-bj.x, gradient);
        pi.dii_boundary+=(-dt*dt*bj.psi/pow(pi.rho,2))*gradient;
    }
}

void System::computeAii( int i)
{
    Particle& pi=particles[i]; pi.aii=0.0;
    Vec2d gradient_ij(0.0,0.0), dji(0.0,0.0);

    for(int& j : pi.fluidNeighbor)
    {
        if(i!=j)
        {
            Particle& pj=particles[j];
            dji=computeDij(j,i);
            p_kernel.monaghanGradient(pi.x-pj.x, gradient_ij);
            pi.aii+=mass*((pi.dii_fluid+pi.dii_boundary)-dji).dot(gradient_ij);
        }
    }
    for(int& j : pi.boundaryNeighbor)
    {
        Boundary& bj=boundaries[j];
        gradient_ij(0.0,0.0);
        p_kernel.monaghanGradient(pi.x-bj.x, gradient_ij);
        pi.aii+=bj.psi*(pi.dii_fluid+pi.dii_boundary).dot(gradient_ij);
    }
}

void System::computeBoundaryVolume()
{
#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(size_t i=0; i<boundaries.size();++i)
    {
        double densityNumber=0.0;
        vector<int> neighbors;
        getNearestNeighbor(neighbors, boundaryGrid, boundaries[i].x);
        for(int& j : neighbors)
        {
            densityNumber += p_kernel.monaghanValue(boundaries[i].x-boundaries[j].x);
        }
        boundaries[i].psi = restDensity/densityNumber;
    }
}

void System::computeMeanDensity()
{
    mean_density=0.0;
    for(int i=0; i<particleNumber; ++i)
    {
        mean_density+=particles[i].rho;
    }
    mean_density/=particleNumber;
}

void System::computeDensityFluctuation()
{
    density_fluctuation=mean_density-restDensity;
}

void System::computeVolume()
{
    real_volume=0.0;
    for(int i=0; i<particleNumber; ++i)
    {
        real_volume += mass/particles[i].rho;
    }
}

void System::addMovingBoundaryBox(Vec2d& offset, Vec2d& scale, Vec2d &translation, bool& oscillatory, Vec2d& delta)
{
    std::vector<int> id = addBoundaryBox(offset, scale);
    m_moving_boundaries.push_back(MovingBoundary(id, translation, oscillatory, delta));
}

std::vector<int> System::addBoundaryBox(Vec2d& offset, Vec2d& scale)
{
    std::vector<int> boundaryId;

    double spacing = h;
    int widthSize = std::floor(scale[0]/spacing);
    int heightSize = std::floor(scale[1]/spacing);
    //Bottom+Top
    for(int i=0; i<=widthSize; ++i)
    {
        Vec2d bPosition(offset[0]+i*spacing, offset[1]);
        boundaries.push_back(Boundary(bPosition,Vec2d(0.0,0.0),0.0));
        boundaryId.push_back(boundaries.size()-1);
        boundaryNumber++;
        Vec2d tPosition(offset[0]+i*spacing, offset[1]+scale[1]);
        boundaries.push_back(Boundary(tPosition,Vec2d(0.0,0.0),0.0));
        boundaryId.push_back(boundaries.size()-1);
        boundaryNumber++;

    }
    //Left+Right
    for(int i=0; i<=heightSize; ++i)
    {
        Vec2d lPosition(offset[0],offset[1]+i*spacing);
        boundaries.push_back(Boundary(lPosition,Vec2d(0.0,0.0),0.0));
        boundaryId.push_back(boundaries.size()-1);
        boundaryNumber++;
        Vec2d rPosition(offset[0]+scale[0],offset[1]+i*spacing);
        boundaries.push_back(Boundary(rPosition,Vec2d(0.0,0.0),0.0));
        boundaryId.push_back(boundaries.size()-1);
        boundaryNumber++;
    }

    double eOffset = 2.0*h;
    double eScale = 4.0*h;
    Vec2d gOffset = offset-Vec2d(eOffset,eOffset);
    Vec2d gScale = scale+Vec2d(eScale,eScale);
    gridInfo.update(gOffset, gScale, 2.0*h);

    return boundaryId;
}

void System::cleanFluidParticle()
{
    particles.clear();
    particleNumber=0;
}

void System::addParticleBox(Vec2d offset, Vec2d scale)
{
    int widthSize = std::floor(scale[0]/h);
    int heightSize = std::floor(scale[1]/h);

    for(int i=0; i<widthSize; ++i)
    {
        for(int j=0; j<heightSize; ++j)
        {
            Vec2d _x = offset + Vec2d(i*h,j*h);
            Vec2d _v(0.0,0.0);
            Vec2d _c(0.0,0.0);
            particles.push_back( Particle(_x,_v,_c) );
            particleNumber++;
        }
    }
}

void System::setGravity(const Vec2d& _gravity)
{
    gravity = _gravity;
}

const Vec2d& System::getGravity() const
{
    return gravity;
}

void System::updateParameters()
{
    ///TODO
}

void System::setParameters( int _wishedNumber, double _volume )
{
    time = 0.0;
    countTime = 0.0;
    countExport = 0;
    particleNumber = 0;
    boundaryNumber = 0;
    volume = _volume;

    maxEta=1.0;
    restDensity = 100;
    _wishedNumber = _wishedNumber*6.5;
    mass = (restDensity * volume) / (_wishedNumber);
    //mass = (restDensity * volume) / (_wishedNumber);
    particlePerCell = 38;
    //h = pow( volume / _wishedNumber, 1.0/2.0);
    //h = 0.5*pow( double(volume*particlePerCell) / double(M_PI*_wishedNumber), 1.0/2.0);
    h = 0.5*pow( double(volume*particlePerCell) / double(M_PI*_wishedNumber), 1.0/2.0);
    double eta = 0.01;
    double H = 0.1;
    double vf = sqrt( 2*9.81*H );
    cs = vf/(sqrt(eta));

    alpha = 1e-3; ///Viscosity
    fcohesion = 1e-3; ///Fluid cohesion, surface tension
    badhesion = 1e-3; ///Boundary adhesion
    sigma=1.0; ///Boundary friction
    boundaryH = h/2.0; //BoundaryH must be <= h (neighbor search purpose)

    dt = 1e-10; ///Time step

    p_kernel = MonaghanKernel( h );
    a_kernel = AkinciKernel( 2.0*h );
    b_kernel = BoundaryKernel( boundaryH, cs );
}

void System::set_benchmark_sink(double timeStep)
{
    this->setTimeStep(timeStep);
    this->setFluidCohesion(0.05);
    this->setViscosity(0.000);
    this->setBoundaryAdhesion(0);
    this->setBoundaryFriction(0);

    Vec2d securityOffset(1.05*this->getSmoothingRadius(),1.05*this->getSmoothingRadius());
    Vec2d boundaryOffset(0,0);
    boundaryOffset -= 1.0*securityOffset;
    Vec2d boundaryBox(1.0,1.0);
    boundaryBox += securityOffset;
    this->addBoundaryBox(boundaryOffset, boundaryBox);

    double startTime, endTime, delay, spacing;
    Vec2d position;
    startTime=0;
    endTime = 4.0;
    delay=0.015;
    spacing = 1.05*this->getSmoothingRadius();
    double velocity = 1.5;
    double scale = 0.2;
    double orientation = M_PI/2.0;
    position = Vec2d(0.0,0.2);
    ParticleSource source1(startTime, endTime, delay, spacing, position, orientation, scale, velocity);
    this->addParticleSource(source1);

    Vec2d minBB(0.95,-0.05);
    Vec2d maxBB(1.1,1.1);
    ParticleSink sink1(minBB, maxBB);
    this->addParticleSink(sink1);
}

void System::set_benchmark_source(double timeStep)
{
    this->setTimeStep(timeStep);
    this->setFluidCohesion(0.05);
    this->setViscosity(0.000);
    this->setBoundaryAdhesion(0);
    this->setBoundaryFriction(0);

    Vec2d securityOffset(1.05*this->getSmoothingRadius(),1.05*this->getSmoothingRadius());
    Vec2d boundaryOffset(0,0);
    boundaryOffset -= 1.0*securityOffset;
    Vec2d boundaryBox(1.0,1.0);
    boundaryBox += securityOffset;
    this->addBoundaryBox(boundaryOffset, boundaryBox);

    double startTime, endTime, delay, spacing;
    Vec2d position;
    startTime=0;
    endTime = 2.0;
    delay=0.02;
    spacing = 1.05*this->getSmoothingRadius();
    double velocity = 1.5;
    double scale = 0.2;
    double orientation = M_PI/2.0;
    position = securityOffset + Vec2d(0.0,0.7);
    ParticleSource source1(startTime, endTime, delay, spacing, position, orientation, scale, velocity);
    this->addParticleSource(source1);
}

void System::set_benchmark_dambreak(double timeStep)
{
    this->setTimeStep(timeStep);
    this->setFluidCohesion(0.05);
    this->setViscosity(0.000);
    this->setBoundaryAdhesion(0);
    this->setBoundaryFriction(0);

    Vec2d securityOffset(1.05*this->getSmoothingRadius(),1.05*this->getSmoothingRadius());
    Vec2d fluidOffset(0.0,0);
    Vec2d fluidBox(0.6,0.6);
    this->addParticleBox(fluidOffset, fluidBox);

    Vec2d boundaryOffset(0,0);
    boundaryOffset -= 1.0*securityOffset;
    Vec2d boundaryBox(1.0,1.0);
    boundaryBox += securityOffset;
    this->addBoundaryBox(boundaryOffset, boundaryBox);
}

void System::set_benchmark_wave(double timeStep)
{
    this->setTimeStep(timeStep);
    this->setFluidCohesion(0);
    this->setViscosity(0);
    this->setBoundaryAdhesion(0);
    this->setBoundaryFriction(0);

    Vec2d securityOffset(1.05*this->getSmoothingRadius(),1.05*this->getSmoothingRadius());
    Vec2d fluidOffset(0.0,0);
    Vec2d fluidBox(1.0,0.2);
    this->addParticleBox(fluidOffset, fluidBox);

    Vec2d boundaryOffset(0,0);
    boundaryOffset -= 1.0*securityOffset;
    Vec2d boundaryBox(1.0,1.0);
    boundaryBox += securityOffset;
    this->addBoundaryBox(boundaryOffset, boundaryBox);

    ///Let the particles place themselves
    init();
    for(size_t i=0; i<100; ++i)
        simulate();

    ///Set moving boundary
    Vec2d movingBoundaryOffset(-1.0,0.0);
    movingBoundaryOffset -= securityOffset;
    Vec2d movingBoundaryScale(1.0,1.0);
    bool oscillatory = true;
    Vec2d delta(0.1,0.1);
    Vec2d boundaryVelocity(0.2,0.0);
    Vec2d translationMotion = this->getTimeStep()*boundaryVelocity;
    this->addMovingBoundaryBox(movingBoundaryOffset, movingBoundaryScale, translationMotion, oscillatory, delta);
}

void System::set_benchmark_splash(double timeStep)
{
    this->setTimeStep(timeStep);
    this->setFluidCohesion(0.1);
    this->setViscosity(0);
    this->setBoundaryAdhesion(0);
    this->setBoundaryFriction(0);

    Vec2d securityOffset(1.05*this->getSmoothingRadius(),1.05*this->getSmoothingRadius());
    Vec2d fluidOffset(0.0,0);
    Vec2d fluidBox(1.0,0.2);
    this->addParticleBox(fluidOffset, fluidBox);

    Vec2d boundaryOffset(0,0);
    boundaryOffset -= 1.0*securityOffset;
    Vec2d boundaryBox(1.0,1.0);
    boundaryBox += securityOffset;
    this->addBoundaryBox(boundaryOffset, boundaryBox);

    ///Let the particles place themselves
    init();
    for(size_t i=0; i<100; ++i)
        simulate();

    ///Add a falling chunk of water
    fluidOffset = Vec2d(0.4,0.3);
    fluidBox = Vec2d(0.2,0.2);
    this->addParticleBox(fluidOffset, fluidBox);
}

void System::init()
{
    prepareGrid();
    computeBoundaryVolume();
    //debugFluid();
}

void System::predictAdvection()
{
#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<particleNumber; ++i)
        computeRho(i);

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<particleNumber; ++i)
        computeNormal(i);

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<particleNumber; ++i)
    {
        computeAdvectionForces(i);
        predictVelocity(i);
        computeDii(i);
    }

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<particleNumber; ++i)
    {
        predictRho(i);
        initializePressure(i);
        computeAii(i);
    }
}

void System::pressureSolve()
{
    int l=0; rho_avg_l = 0.0;

    while( ( (rho_avg_l-restDensity)>maxEta ) || (l<2) )
    {
#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
        for(int i=0; i<particleNumber; ++i)
            computeSumDijPj(i);

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
        for(int i=0; i<particleNumber; ++i)
        {
            computePressure(i);
        }

        computeError();

        ++l;

        //debugIteration(l);
    }

}

void System::integration()
{
    countTime++; time+=dt;

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<particleNumber; ++i)
        computePressureForce(i);

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<particleNumber; ++i)
    {
        Particle& pi=particles[i];
        pi.v = pi.v_adv + (dt*pi.f_p)/mass;
        pi.x += dt*pi.v;
    }
}

void System::prepareGrid()
{
    updateGridInfo();

    //Update fluid grid
    fluidGrid.resize(gridInfo.size());
    std::vector<int> defaultVec; defaultVec.reserve(60);
    std::fill(fluidGrid.begin(), fluidGrid.end(), defaultVec);

    for(size_t i=0; i<particles.size(); ++i)
    {
        int id = gridInfo.cellId(particles[i].x);
        if(gridInfo.isInside(id))
            fluidGrid[id].push_back(i);
    }

    //Update boundary grid
    boundaryGrid.resize(gridInfo.size());
    std::fill(boundaryGrid.begin(), boundaryGrid.end(), std::vector<int>());

    for(size_t i=0; i<boundaries.size(); ++i)
    {
        int id = gridInfo.cellId(boundaries[i].x);
        if(gridInfo.isInside(id))
            boundaryGrid[id].push_back(i);
    }

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(size_t i = 0; i < particles.size(); ++i)
        getNearestNeighbor(i, 2.0*h);
}

void System::getNearestNeighbor(const int i, const double radius)
{
    Particle& p = particles[i];
    p.fluidNeighbor.clear();
    p.fluidNeighbor.reserve(60);
    p.boundaryNeighbor.clear();
    p.boundaryNeighbor.reserve(60);

    std::vector<int> neighborCell;
    gridInfo.get9Neighbors(neighborCell, p.x, radius);

    Vec2d d(0.0,0.0);
    int bParticleId=-1;
    int fParticleId=-1;
    for(size_t i=0; i<neighborCell.size(); ++i)
    {
        std::vector<int>& bNeighborCell = boundaryGrid[neighborCell[i]];
        for(size_t j=0; j<bNeighborCell.size(); ++j)
        {
            bParticleId = boundaryGrid[neighborCell[i]][j];
            Boundary& bParticle = boundaries[bParticleId];
            d = bParticle.x-p.x;
            if( d.squaredNorm()<radius*radius )
                p.boundaryNeighbor.push_back(bParticleId);
        }

        std::vector<int>& fNeighborCell = fluidGrid[neighborCell[i]];
        for(size_t j=0; j< fNeighborCell.size(); ++j)
        {
            fParticleId = fluidGrid[neighborCell[i]][j];
            Particle& fParticle = particles[fParticleId];
            d = fParticle.x-p.x;
            if( d.squaredNorm()<radius*radius)
                p.fluidNeighbor.push_back(fParticleId);
        }
    }
}

void System::getNearestNeighbor(vector< int >& neighbor, const vector< vector<int> >& grid, const Vec2d &x)
{
    std::vector<int> neighborCell;
    gridInfo.get9Neighbors(neighborCell, x, gridInfo.spacing());
    neighbor.clear();
    for(size_t i=0; i<neighborCell.size(); ++i)
    {
        for(size_t j=0; j<grid[neighborCell[i]].size(); ++j)
        {
            neighbor.push_back(grid[neighborCell[i]][j]);
        }
    }
}

void System::getNearestFluidNeighbor(vector< int >& neighbor, const Vec2d &x, const double radius)
{
    neighbor.clear();

    Vec2d d(0.0,0.0);
    int fParticleId=-1;
    std::vector<int> neighborCell;
    gridInfo.get9Neighbors(neighborCell, x, gridInfo.spacing());

    for(size_t i=0; i<neighborCell.size(); ++i)
    {
        std::vector<int>& fNeighborCell = fluidGrid[neighborCell[i]];
        for(size_t j=0; j<fNeighborCell.size(); ++j)
        {
            fParticleId = fluidGrid[neighborCell[i]][j];
            Particle& fParticle = particles[fParticleId];
            d = fParticle.x-x;
            if( d.squaredNorm()<radius*radius)
                neighbor.push_back(fParticleId);
        }
    }
}

void System::simulate()
{
    prepareGrid();
    computeBoundaryVolume();
    predictAdvection();
    pressureSolve();
    integration();
    applyMovingBoundaries();
    applySources();
    applySinks();
    computeStats();
}

void System::applyMovingBoundaries()
{
    for(MovingBoundary& movingBoundary : m_moving_boundaries)
    {
        movingBoundary.applyMotion(boundaries);
    }
}

void System::computeStats()
{
    computeMeanDensity();
    computeVolume();
    computeDensityFluctuation();
}

void System::debugIteration(int l)
{
    std::cout.precision(10);
    std::cout << "rest density " << restDensity << std::endl;
    std::cout << "rho avg : " << rho_avg_l << std::endl;
    std::cout << "l : " << l << std::endl;
}

void System::debugFluid()
{
    std::cout << "Particle Number : " << particleNumber << std::endl;
    std::cout << "Boundary Number : " << boundaryNumber << std::endl;
    std::cout << "Smoothing Radius : " << h << std::endl;
    std::cout << "Radius : " << pow( 3.0*volume/(4*M_PI*particleNumber), 1.0/3.0 ) << std::endl;
    std::cout << "Speed sound : " << cs << std::endl;
    std::cout << "Timestep : " << dt << std::endl;
    std::cout << std::endl;
}

}
