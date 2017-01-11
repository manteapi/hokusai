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
#include "../include/hokusai/rasterizer.hpp"
#include "../include/hokusai/gridUtility.hpp"
#include "../include/hokusai/HSL2RGB.hpp"
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
    m_fluidParams = FluidParams();
    m_boundaryParams = BoundaryParams();
    m_solverParams = SolverParams();
    m_countExport = 0;
    m_countTime = 0;
    m_particleNumber = 0;
    m_boundaryNumber = 0;

    m_meanDensity = 0.0;
    m_densityFluctuation = 0.0;
    m_realVolume = 0.0;
    m_time = 0.0;
    m_gravity = Vec3r(0,-9.81,0);
    m_gridInfo = GridUtility();
    m_particles = std::vector<Particle>();
    m_boundaries = std::vector<Boundary>();
}

System::System(const FluidParams &fluidParams, const BoundaryParams& boundaryParams, const SolverParams& solverParams)
{
    m_fluidParams = fluidParams;
    m_boundaryParams = boundaryParams;
    m_solverParams = solverParams;
    m_countExport = 0;
    m_countTime = 0;
    m_particleNumber = 0;
    m_boundaryNumber = 0;

    m_meanDensity = 0.0;
    m_densityFluctuation = 0.0;
    m_realVolume = 0.0;
    m_time = 0.0;

    m_gravity = Vec3r(0,-9.81,0);

    m_gridInfo = GridUtility();
    m_particles = std::vector<Particle>();
    m_boundaries = std::vector<Boundary>();
}

System::~System(){}

void System::computeDensity(int i)
{
    Particle& pi=m_particles[i];
    std::vector<int>& fneighbors=pi.fluidNeighbor;
    std::vector<int>& bneighbors=pi.boundaryNeighbor;
    pi.rho=0.0;
    for(int& j : fneighbors)
    {
        Particle& pj=m_particles[j];
        pi.rho += m_fluidParams.mass()*m_fluidParams.monaghanKernel().monaghanValue(pi.x-pj.x);
    }
    for(int& j : bneighbors)
    {
        Boundary& bj = m_boundaries[j];
        pi.rho += m_fluidParams.monaghanKernel().monaghanValue(pi.x-bj.x)*bj.psi;
    }
}

void System::computeNormal(int i)
{
    //Compute normal
    Particle& pi=m_particles[i];
    std::vector< int > & neighbors = m_particles[i].fluidNeighbor;
    Vec3r n(0.0);
    for(size_t j=0; j<neighbors.size(); ++j)
    {
        if(i!=neighbors[j])
        {
            Particle& pj = m_particles[neighbors[j]];
            n += (m_fluidParams.mass()/pj.rho)*m_fluidFluidMonaghanGradient[i][j];
        }
    }
    pi.n = m_fluidParams.smoothingRadius()*n;
}

bool System::isSurfaceParticle(int i, HReal treshold)
{
    HReal n_length = m_particles[i].n.lengthSquared();
    if( n_length > treshold )
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::vector<Particle> System::getSurfaceParticle()
{
#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<m_particleNumber; ++i)
        computeDensity(i);

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<m_particleNumber; ++i)
        computeNormal(i);

    HReal treshold = 0.05;
    std::vector<Particle> surfaceParticles;
    for(int i=0; i<m_particleNumber; ++i)
    {
        if( isSurfaceParticle(i, treshold) )
            surfaceParticles.push_back(m_particles[i]);
    }
    return surfaceParticles;
}

void System::predictVelocity(int i)
{
    Particle& pi=m_particles[i];
    pi.v_adv = pi.v + (m_solverParams.timeStep()/m_fluidParams.mass())*pi.f_adv;
}

void System::predictDensity(int i)
{
    Particle& pi=m_particles[i];
    std::vector<int>& fneighbors=pi.fluidNeighbor;
    std::vector<int>& bneighbors=pi.boundaryNeighbor;
    HReal fdrho=0.0, bdrho=0.0;
    Vec3r vij_adv(0.0);

    for(size_t j=0; j<fneighbors.size(); ++j)
    {
        if(i!=fneighbors[j])
        {
            Particle& pj=m_particles[fneighbors[j]];
            vij_adv=pi.v_adv-pj.v_adv;
            fdrho+=m_fluidParams.mass()*Vec3r::dotProduct(vij_adv, m_fluidFluidMonaghanGradient[i][j]);
        }
    }

    Vec3r vb(0.1), v(0.0);
    for(size_t j=0; j<bneighbors.size(); ++j)
    {
        Boundary& bj=m_boundaries[bneighbors[j]];
        v = pi.v_adv - vb; //vb(t+dt)
        bdrho+=bj.psi*Vec3r::dotProduct(v,m_fluidBoundaryMonaghanGradient[i][j]);
    }

    pi.rho_adv = pi.rho + m_solverParams.timeStep()*( fdrho + bdrho );
}

void System::computeSumDijPj(int i)
{
    Particle& pi=m_particles[i];
    std::vector<int>& fneighbors=pi.fluidNeighbor;
    pi.sum_dij.fill(0.0);
    for(size_t j=0; j<fneighbors.size(); ++j)
    {
        if(i!=fneighbors[j])
        {
            Particle& pj=m_particles[fneighbors[j]];
            pi.sum_dij+=(-m_fluidParams.mass()/(pj.rho*pj.rho))*pj.p_l*m_fluidFluidMonaghanGradient[i][j];
        }
    }
    pi.sum_dij *= m_solverParams.timeStep()*m_solverParams.timeStep();
}

void System::computeViscosityForces(const int& i)
{
    HReal kij=0.0, epsilon=0.01, Pij=0.01, dotVijRij=0.0;
    Vec3r r(0.0), vij(0.0);
    Particle& pi=m_particles[i];
    for(size_t j=0; j<pi.fluidNeighbor.size(); ++j)
    {
        Particle& pj=m_particles[pi.fluidNeighbor[j]];
        r = pi.x - pj.x;
        vij = pi.v - pj.v;
        dotVijRij = Vec3r::dotProduct(vij,r);
        if(dotVijRij < 0)
        {
            kij = (m_fluidParams.restDensity()+m_fluidParams.restDensity())/(pi.rho+pj.rho);
            epsilon=0.01;
            Pij = -kij*(2.0*m_fluidParams.viscosity()*m_fluidParams.smoothingRadius()*m_fluidParams.soundSpeed()/(pi.rho+pj.rho)) * ( dotVijRij / (r.lengthSquared() + epsilon*m_fluidParams.smoothingRadius()*m_fluidParams.smoothingRadius()) );
            pi.f_adv += -kij*m_fluidParams.mass()*m_fluidParams.mass()*Pij*m_fluidFluidMonaghanGradient[i][j];
        }
    }
}

void System::computeSurfaceTensionForces(const int & i)
{
    HReal kij=0.0, l=0.0;
    Vec3r r(0.0), nij(0.0), cohesionForce(0.0), curvatureForce(0.0);
    Particle& pi=m_particles[i];
    for(size_t j=0; j<pi.fluidNeighbor.size(); ++j)
    {
        if(i!=pi.fluidNeighbor[j])
        {
            Particle& pj=m_particles[pi.fluidNeighbor[j]];
            if(pi.isSurface==true || pj.isSurface==true)
            {
                r = pi.x - pj.x;
                kij = (m_fluidParams.restDensity()+m_fluidParams.restDensity())/(pi.rho+pj.rho);
                l = r.length();
                cohesionForce = -(m_fluidParams.cohesion()*m_fluidParams.mass()*m_fluidParams.mass()*m_fluidParams.akinciKernel().cohesionValue(l)/l) * r;
                nij = pi.n-pj.n;
                curvatureForce = -m_fluidParams.cohesion()*m_fluidParams.mass()*nij;
                pi.f_adv += kij*(cohesionForce+curvatureForce);
            }
        }
    }
}


void System::computeBoundaryFrictionForces(const int& i)
{
    HReal epsilon=0.01, nu=0.0, Pij=0.0, dotVijRij=0.0;
    Vec3r vij(0.0), xij(0.0);
    Particle& pi=m_particles[i];
    for(size_t j=0; j<pi.boundaryNeighbor.size(); ++j)
    {
        Boundary& bj=m_boundaries[pi.boundaryNeighbor[j]];
        vij = pi.v;//-pj.v;
        xij= pi.x - bj.x;
        dotVijRij = Vec3r::dotProduct(vij,xij);
        if(dotVijRij<0)
        {
            epsilon=0.01;
            nu = (m_boundaryParams.friction()*m_fluidParams.smoothingRadius()*m_fluidParams.soundSpeed())/(2.0*pi.rho);
            Pij = -nu * ( std::min(dotVijRij, (HReal)(0.0)) / (xij.lengthSquared() + epsilon*m_fluidParams.smoothingRadius()*m_fluidParams.smoothingRadius()) );
            pi.f_adv += -m_fluidParams.mass()*bj.psi*Pij*m_fluidBoundaryMonaghanGradient[i][j];
        }
    }
}

void System::computeBoundaryAdhesionForces(const int& i)
{
    HReal l=0.0;
    Vec3r xij(0.0);
    Particle& pi=m_particles[i];
    for(size_t j=0; j<pi.boundaryNeighbor.size(); ++j)
    {
        Boundary& bj=m_boundaries[pi.boundaryNeighbor[j]];
        xij= pi.x - bj.x;
        l = xij.length();
        pi.f_adv += -(m_boundaryParams.adhesion()*m_fluidParams.mass()*m_boundaries[j].psi*m_fluidParams.akinciKernel().adhesionValue(l)/l)*xij;
    }
}

void System::computePressure(int i)
{
    Particle& pi=m_particles[i];
    std::vector<int>& fneighbors=pi.fluidNeighbor, bneighbors=pi.boundaryNeighbor;
    HReal fsum=0.0, bsum=0.0, omega=0.5;
    Vec3r dji(0.0), aux(0.0), r(0.0);

    for(size_t j=0; j<fneighbors.size(); ++j)
    {
        if(i!=fneighbors[j])
        {
            Particle& pj=m_particles[fneighbors[j]];
            //Compute dji
            dji=-(m_solverParams.timeStep()*m_solverParams.timeStep()*m_fluidParams.mass())/(pi.rho*pi.rho)*(-m_fluidFluidMonaghanGradient[i][j]);
            //Compute fsum
            aux = pi.sum_dij - (pj.dii_fluid+pj.dii_boundary)*pj.p_l - (pj.sum_dij - dji*pi.p_l);
            fsum+=m_fluidParams.mass()*Vec3r::dotProduct(aux, m_fluidFluidMonaghanGradient[i][j]);
        }
    }

    for(size_t j=0; j<bneighbors.size(); ++j)
    {
        Boundary& bj=m_boundaries[bneighbors[j]];
        r=pi.x-bj.x;
        bsum+=bj.psi*Vec3r::dotProduct(pi.sum_dij,m_fluidBoundaryMonaghanGradient[i][j]);
    }

    HReal previousPl = pi.p_l;
    pi.rho_corr = pi.rho_adv + fsum + bsum;
    if(std::abs(pi.aii)>std::numeric_limits<HReal>::epsilon())
    {
        pi.p_l = (1-omega)*previousPl + (omega/pi.aii)*(m_fluidParams.restDensity() - pi.rho_corr);
    }
    else
    {
        pi.p_l = 0.0;
    }
    pi.p = std::max(pi.p_l, (HReal)(0.0));
    pi.p_l = pi.p;
    pi.rho_corr += pi.aii*previousPl;
}

void System::computePressureForce(int i)
{
    Particle& pi=m_particles[i];
    pi.f_p.fill(0.0);
    std::vector<int>& fneighbors=pi.fluidNeighbor;
    std::vector<int>& bneighbors=pi.boundaryNeighbor;

    //Fluid Pressure Force
    for(size_t j=0; j<fneighbors.size(); ++j)
    {
        computeFluidPressureForce(i, fneighbors[j], m_fluidFluidMonaghanGradient[i][j]);
    }

    //Boundary Pressure Force [Akinci 2012]
    for(size_t j=0; j<bneighbors.size(); ++j)
    {
        computeBoundaryPressureForce(i, bneighbors[j], m_fluidBoundaryMonaghanGradient[i][j]);
    }
}

void System::computeFluidPressureForce(int i, int j, const Vec3r& gradient_ij)
{
    Particle& pi=m_particles[i];
    Particle& pj=m_particles[j];
    if( i!=j )
    {
        pi.f_p += -m_fluidParams.mass()*m_fluidParams.mass()*( pi.p/(pi.rho*pi.rho) + pj.p/(pj.rho*pj.rho) ) * gradient_ij;
    }
}

void System::computeBoundaryPressureForce(int i, int j, const Vec3r& gradient_ij)
{
    Particle& pi=m_particles[i];
    Boundary& bj=m_boundaries[j];
    pi.f_p += -m_fluidParams.mass()*bj.psi*( pi.p/(pi.rho*pi.rho) ) * gradient_ij;
}

void System::initializePressure(int i)
{
    Particle& pi=m_particles[i];
    pi.p_l=0.5*pi.p;
}

void System::computeError()
{
    m_solverParams.averageDensity() = 0.0;
    for(int i=0; i<m_particleNumber; ++i)
        m_solverParams.averageDensity()+=m_particles[i].rho_corr;
    m_solverParams.averageDensity() /= m_particleNumber;
}

void System::computeDii(int i)
{
    Particle& pi=m_particles[i];
    pi.dii_fluid.fill(0.0);
    pi.dii_boundary.fill(0.0);
    for(size_t j=0; j<pi.fluidNeighbor.size(); ++j)
    {
        if(i!=pi.fluidNeighbor[j])
        {
            Particle& pj=m_particles[pi.fluidNeighbor[j]];
            pi.dii_fluid+=(-m_solverParams.timeStep()*m_solverParams.timeStep()*m_fluidParams.mass()/(pi.rho*pi.rho))*m_fluidFluidMonaghanGradient[i][j];
        }
    }
    for(size_t j=0; j<pi.boundaryNeighbor.size(); ++j)
    {
        Boundary& bj=m_boundaries[pi.boundaryNeighbor[j]];
        pi.dii_boundary+=(-m_solverParams.timeStep()*m_solverParams.timeStep()*bj.psi/(pi.rho*pi.rho))*m_fluidBoundaryMonaghanGradient[i][j];
    }
}

void System::computeAii( int i)
{
    Particle& pi=m_particles[i]; pi.aii=0.0;
    Vec3r dji(0.0);
    for(size_t j=0; j<pi.fluidNeighbor.size(); ++j)
    {
        if(i!=pi.fluidNeighbor[j])
        {
            Particle& pj=m_particles[pi.fluidNeighbor[j]];
            //Compute dji
            dji=-(m_solverParams.timeStep()*m_solverParams.timeStep()*m_fluidParams.mass())/(pi.rho*pi.rho)*(-m_fluidFluidMonaghanGradient[i][j]);
            //Compute aii
            pi.aii+=m_fluidParams.mass()*Vec3r::dotProduct((pi.dii_fluid+pi.dii_boundary)-dji,m_fluidFluidMonaghanGradient[i][j]);
        }
    }
    for(size_t j=0; j<pi.boundaryNeighbor.size(); ++j)
    {
        Boundary& bj=m_boundaries[pi.boundaryNeighbor[j]];
        pi.aii+=bj.psi*Vec3r::dotProduct(pi.dii_fluid+pi.dii_boundary,m_fluidBoundaryMonaghanGradient[i][j]);
    }
}

void System::getNearestNeighbor(const int& particleId, const HReal& radius)
{
    Particle& p = m_particles[particleId];

    size_t lastFluidSize=p.fluidNeighbor.size();
    p.fluidNeighbor.clear();
    p.fluidNeighbor.reserve(lastFluidSize);

    size_t lastBoundarySize=p.boundaryNeighbor.size();
    p.boundaryNeighbor.clear();
    p.boundaryNeighbor.reserve(lastBoundarySize);

    HReal radiusSquared = radius*radius;
    int bParticleId = 0, fParticleId=0;
    Vec3r d(0.0);

    std::vector<int> neighborCell;
    m_gridInfo.get27Neighbors(neighborCell, p.x, radius);

    for(size_t i=0; i<neighborCell.size(); ++i)
    {
        std::vector<int>& bNeighborCell = m_boundaryGrid[neighborCell[i]];
        for(size_t j=0; j<bNeighborCell.size(); ++j)
        {
            bParticleId = m_boundaryGrid[neighborCell[i]][j];
            Boundary& bParticle = m_boundaries[bParticleId];
            d = bParticle.x-p.x;
            if(d.lengthSquared()<radiusSquared)
            {
                p.boundaryNeighbor.push_back(bParticleId);
            }
        }

        std::vector<int>& fNeighborCell = m_fluidGrid[neighborCell[i]];
        for(size_t j=0; j< fNeighborCell.size(); ++j)
        {
            fParticleId = m_fluidGrid[neighborCell[i]][j];
            Particle& fParticle = m_particles[fParticleId];
            d = fParticle.x-p.x;
            if(d.lengthSquared()<radiusSquared)
            {
                p.fluidNeighbor.push_back(fParticleId);
            }
        }
    }
}

BoundaryParams& System::boundaryParams()
{
    return m_boundaryParams;
}

const BoundaryParams& System::boundaryParams() const
{
    return m_boundaryParams;
}

FluidParams& System::fluidParams()
{
    return m_fluidParams;
}

const FluidParams& System::fluidParams() const
{
    return m_fluidParams;
}

SolverParams& System::solverParams()
{
    return m_solverParams;
}

const SolverParams& System::solverParams() const
{
    return m_solverParams;
}

std::vector<Particle>& System::particles()
{
    return m_particles;
}

const std::vector<Particle>& System::particles() const
{
    return m_particles;
}

std::vector<Boundary>& System::boundaries()
{
    return m_boundaries;
}

const std::vector<Boundary>& System::boundaries() const
{
    return m_boundaries;
}

void System::getNearestNeighbor(std::vector< int >& neighbor, const std::vector< std::vector<int> >& grid, const Vec3r &x)
{
    std::vector<int> neighborCell;
    m_gridInfo.get27Neighbors(neighborCell, x, m_gridInfo.spacing());

    size_t lastSize = neighbor.size();
    neighbor.clear();
    neighbor.reserve(lastSize);
    for(size_t i=0; i<neighborCell.size(); ++i)
    {
        for(size_t j=0; j<grid[neighborCell[i]].size(); ++j)
        {
            neighbor.push_back(grid[neighborCell[i]][j]);
        }
    }
}

void System::computeBoundaryVolume()
{
    for(int i=0; i<m_boundaryNumber; ++i)
    {
        HReal densityNumber=0.0;
        std::vector<int> neighbors;
        getNearestNeighbor(neighbors, m_boundaryGrid, m_boundaries[i].x);
        for(int& j : neighbors)
            densityNumber += m_fluidParams.monaghanKernel().monaghanValue(m_boundaries[i].x-m_boundaries[j].x);
        m_boundaries[i].psi = m_fluidParams.restDensity()/densityNumber;
    }
}

void System::computeMeanDensity()
{
    m_meanDensity=0.0;
    for(int i=0; i<m_particleNumber; ++i)
    {
        m_meanDensity+=m_particles[i].rho;
    }
    m_meanDensity/=m_particleNumber;
}

void System::computeDensityFluctuation()
{
    m_densityFluctuation=m_meanDensity-m_fluidParams.restDensity();
}

void System::computeVolume()
{
    m_realVolume=0.0;
    for(int i=0; i<m_particleNumber; ++i)
    {
        Particle& pi = m_particles[i];
        m_realVolume += m_fluidParams.mass()/pi.rho;
    }
}

void System::setGravity(const Vec3r& _gravity)
{
    m_gravity = _gravity;
}

const Vec3r& System::getGravity()
{
    return m_gravity;
}

void System::init()
{  
    prepareGrid();

    computeBoundaryVolume();

    for(size_t i=0; i<m_particles.size(); ++i)
    {
        m_particles[i].isSurface = true;
    }

    debugFluid();
}

void System::addBoundaryBox(const Vec3r& offset, const Vec3r& scale)
{    
    std::vector<Vec3r> positions = getBoxSampling(offset, scale, m_fluidParams.smoothingRadius());
    for(const Vec3r& x : positions)
    {
        m_boundaries.push_back(Boundary(x,Vec3r(0.0),0.0));
        m_boundaryNumber++;
    }
}

void System::addBoundarySphere(const Vec3r& offset, const HReal& radius)
{
    std::vector<Vec3r> samples = getSphereSampling(offset, radius, m_fluidParams.smoothingRadius(), m_fluidParams.smoothingRadius());
    for(size_t i=0; i<samples.size(); ++i)
    {
        m_boundaries.push_back(Boundary(samples[i],Vec3r(0.0),0.0));
        m_boundaryNumber++;
    }
}

void System::addBoundaryHemiSphere(const Vec3r& offset, const HReal& radius)
{
    std::vector<Vec3r> samples = getHemiSphereSampling(offset, radius, m_fluidParams.smoothingRadius(), m_fluidParams.smoothingRadius());
    for(size_t i=0; i<samples.size(); ++i)
    {
        m_boundaries.push_back(Boundary(samples[i],Vec3r(0.0),0.0));
        m_boundaryNumber++;
    }
}

void System::addBoundaryDisk(const Vec3r& offset, const HReal& radius)
{
    std::vector<Vec3r> samples = getDiskSampling(offset, radius, m_fluidParams.smoothingRadius());
    for(size_t i=0; i<samples.size(); ++i)
    {
        m_boundaries.push_back(Boundary(samples[i],Vec3r(0.0),0.0));
        m_boundaryNumber++;
    }
}

void System::addBoundaryCylinder(const Vec3r& offset, const HReal& radius, const HReal& height)
{
    std::vector<Vec3r> samples = getCylinderSampling(offset, height, radius, m_fluidParams.smoothingRadius(), m_fluidParams.smoothingRadius());
    for(size_t i=0; i<samples.size(); ++i)
    {
        m_boundaries.push_back(Boundary(samples[i],Vec3r(0.0),0.0));
        m_boundaryNumber++;
    }
}

void System::translateBoundaries(const Vec3r& t)
{
    for(size_t i=0; i<m_boundaries.size(); ++i)
    {
        m_boundaries[i].x += t;
    }
}

void System::translateParticles(const Vec3r& t)
{
    for(size_t i=0; i<m_particles.size(); ++i)
    {
        m_particles[i].x += t;
    }
}

void System::addParticleSphere(const Vec3r& centre, const HReal radius, const Vec3r& velocity)
{
    std::vector<Vec3r> positions = getBallSampling(centre, radius, m_fluidParams.smoothingRadius());
    for(Vec3r& x : positions)
    {
        m_particles.push_back( Particle(x,velocity) );
        m_particleNumber++;
    }
}

void System::addParticleSource(const ParticleSource& s)
{
    m_pSources.push_back(s);
}

void System::addParticleBox(const Vec3r& offset, const Vec3r& scale, const Vec3r& velocity)
{
    std::vector<Vec3r> positions = getCubeSampling(offset, scale, m_fluidParams.smoothingRadius());
    for(Vec3r & x : positions)
    {
        m_particles.push_back( Particle(x, velocity) );
        m_particleNumber++;
    }
}

void System::addFluidParticle(const Vec3r& x, const Vec3r& v)
{
    m_particles.push_back( Particle(x,v) );
    m_particleNumber++;
}

void System::addBoundaryParticle(const Vec3r& x, const Vec3r& v)
{
    m_boundaries.push_back( Boundary(x,v) );
    m_boundaryNumber++;
}

void System::addBoundaryMesh(const char* filename)
{
    TriMesh mesh(filename);
    std::vector<Vec3r> samples;
    AkinciMeshSampling(mesh, m_fluidParams.smoothingRadius()/2.0, samples);
    for(size_t i=0; i<samples.size(); ++i)
    {
        m_boundaries.push_back(Boundary(samples[i],Vec3r(0.0),0.0));
        m_boundaryNumber++;
    }
}

bool pairCompare( const std::pair<int,int>& e1, const std::pair<int,int>& e2 )
{
    return (e1.second < e2.second);
}

void System::mortonSortBoundary()
{
}

void System::mortonSortFluid()
{
    std::vector< std::pair<int, int> > particleZindex;
    particleZindex.resize(m_particleNumber);
    //Fill particleZindex with particle index and Z-index
    for(int i = 0; i < m_particleNumber; ++i)
    {
        Vec3i _gridIndex = m_gridInfo.worldToGrid(m_particles[i].x);
        std::array<int,3> gridIndex;
        for(int j=0; j<3; ++j)
            gridIndex[j] = _gridIndex[j];
        int zindex = mortonNumber( gridIndex );
        std::pair<int,int> paire(i, zindex);
        particleZindex[i] = paire;
    }

    //Sort according to z-index
    std::sort( particleZindex.begin(), particleZindex.end(), pairCompare );

    //Move particles according to z-index
    std::vector< Particle > oldParticles = m_particles;

    for(int i = 0; i < m_particleNumber; ++i)
    {
        std::pair<int,int>& paire = particleZindex[i];
        if(i != paire.first)
        {
            m_particles[i] = oldParticles[paire.first];
        }
    }
}

void System::computeSurfaceParticle()
{
    for(size_t i=0; i<m_particles.size(); ++i)
    {
        m_particles[i].isSurface = false;
    }

    std::vector<int> tmpParticleStack;
    for(size_t i=0; i<m_particles.size(); ++i)
    {
        //Not good enough
        if( isSurfaceParticle(i, 0.2) || m_particles[i].fluidNeighbor.size() < 0.5*m_particlePerCell)
            tmpParticleStack.push_back(i);
    }

    std::set<int> surfaceParticles;
    for(size_t i=0; i<tmpParticleStack.size(); ++i)
    {
        surfaceParticles.insert( tmpParticleStack[i] );
        for(size_t j=0; j<m_particles[ tmpParticleStack[i] ].fluidNeighbor.size(); ++j)
        {
            surfaceParticles.insert( m_particles[ tmpParticleStack[i] ].fluidNeighbor[j] );
        }
    }

    for(int pId : surfaceParticles)
    {
        m_particles[pId].isSurface=true;
    }

    int surfaceParticle = 0;
    for(size_t i=0; i<m_particles.size(); ++i)
    {
        if(m_particles[i].isSurface == true)
        {
            surfaceParticle++;
        }
    }
}

void System::precomputeKernel()
{
    //Fluid - fluid monaghan gradient
    m_fluidFluidMonaghanGradient.resize( m_particles.size() );
#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(size_t i=0; i<m_particles.size(); ++i)
    {
        const Particle& pi = m_particles[i];
        m_fluidFluidMonaghanGradient[i].resize(pi.fluidNeighbor.size());
        for(size_t j=0; j<pi.fluidNeighbor.size(); ++j)
        {
            const Particle& pj = m_particles[pi.fluidNeighbor[j]];
            if(i!=j)
            {
                m_fluidParams.monaghanKernel().monaghanGradient(pi.x-pj.x, m_fluidFluidMonaghanGradient[i][j]);
            }
            else
            {
                m_fluidFluidMonaghanGradient[i][j] = Vec3r(0.0,0.0,0.0);
            }
        }
    }

    //Fluid - boundary monaghan gradient
    m_fluidBoundaryMonaghanGradient.resize( m_particles.size() );
#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(size_t i=0; i<m_particles.size(); ++i)
    {
        const Particle& pi = m_particles[i];
        m_fluidBoundaryMonaghanGradient[i].resize(pi.boundaryNeighbor.size());
        for(size_t j=0; j<pi.boundaryNeighbor.size(); ++j)
        {
            const Boundary& bj = m_boundaries[pi.boundaryNeighbor[j]];
            m_fluidParams.monaghanKernel().monaghanGradient(pi.x-bj.x, m_fluidBoundaryMonaghanGradient[i][j]);
        }
    }
}

void System::prepareGrid()
{
    if( m_countTime%100 == 0 )
    {
        mortonSortFluid();
        mortonSortBoundary();
    }

    Vec3r minBB(std::numeric_limits<HReal>::max());
    Vec3r maxBB(-std::numeric_limits<HReal>::max());
    for(size_t i=0; i<m_particles.size(); ++i)
    {
        Particle & pi = m_particles[i];
        for(size_t j=0; j<3; ++j)
        {
            minBB[j] = std::min(minBB[j], pi.x[j]);
            maxBB[j] = std::max(maxBB[j], pi.x[j]);
        }
    }
    for(size_t i=0; i<m_boundaries.size(); ++i)
    {
        Boundary & bi = m_boundaries[i];
        for(size_t j=0; j<3; ++j)
        {
            minBB[j] = std::min(minBB[j], bi.x[j]);
            maxBB[j] = std::max(maxBB[j], bi.x[j]);
        }
    }

    Vec3r offset = minBB-Vec3r(2.0*m_fluidParams.smoothingRadius());
    Vec3r scale = maxBB-minBB+Vec3r(4.0*m_fluidParams.smoothingRadius());
    m_gridInfo = GridUtility(offset, scale, 2.0*m_fluidParams.smoothingRadius());

    m_fluidGrid.resize(m_gridInfo.size());
    m_boundaryGrid.resize(m_gridInfo.size());
    for(int i=0; i<m_gridInfo.size(); ++i)
    {
        size_t lastFluidGridSize = m_fluidGrid[i].size();
        m_fluidGrid[i].clear();
        m_fluidGrid[i].reserve(lastFluidGridSize);

        size_t lastBoundaryGridSize = m_boundaryGrid[i].size();
        m_boundaryGrid[i].clear();
        m_boundaryGrid[i].reserve(lastBoundaryGridSize);
    }

    for(size_t i=0; i<m_particles.size(); ++i)
    {
        int id = m_gridInfo.cellId(m_particles[i].x);
        if(m_gridInfo.isInside(id))
            m_fluidGrid[id].push_back(i);
    }

    for(size_t i=0; i<m_boundaries.size(); ++i)
    {
        int id = m_gridInfo.cellId(m_boundaries[i].x);
        if(m_gridInfo.isInside(id))
            m_boundaryGrid[id].push_back(i);
    }

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < m_particleNumber; ++i)
        getNearestNeighbor(i, 2.0*m_fluidParams.smoothingRadius());
}

void System::predictAdvection()
{
#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<m_particleNumber; ++i)
        computeDensity(i);

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<m_particleNumber; ++i)
        computeNormal(i);

    computeSurfaceParticle();

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<m_particleNumber; ++i)
    {
        //Compute advection forces
        m_particles[i].f_adv.fill(0.0);
        computeViscosityForces(i);
        computeSurfaceTensionForces(i);
        computeBoundaryFrictionForces(i);
        computeBoundaryAdhesionForces(i);
        m_particles[i].f_adv+=m_gravity*m_fluidParams.mass();

        predictVelocity(i);

        computeDii(i);
    }

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<m_particleNumber; ++i)
    {
        predictDensity(i);
        initializePressure(i);
        computeAii(i);
    }
}

void System::pressureSolve()
{
    int l=0; m_solverParams.averageDensity() = 0.0;

    while(((m_solverParams.averageDensity()-m_fluidParams.restDensity())>m_solverParams.maxDensityError() ) ||
          (l<m_solverParams.maxPressureSolveIterationNb()) )
    {
#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
        for(int i=0; i<m_particleNumber; ++i)
        {
            computeSumDijPj(i);
        }

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
        for(int i=0; i<m_particleNumber; ++i)
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
    m_countTime++; m_time+=m_solverParams.timeStep();

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<m_particleNumber; ++i)
    {
        computePressureForce(i);
    }

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<m_particleNumber; ++i)
    {
        Particle& pi=m_particles[i];
        pi.v = pi.v_adv + (m_solverParams.timeStep()*pi.f_p)/m_fluidParams.mass();
        pi.x += m_solverParams.timeStep()*pi.v;
    }
}

void System::computeSimulationStep()
{
    prepareGrid();
    precomputeKernel();
    predictAdvection();
    pressureSolve();
    integration();

    applySources();
    applySinks();
    computeStats();
}

void System::applySources()
{
    for(ParticleSource& s : m_pSources)
    {
        std::vector<Particle> p_new = s.apply(this->m_time);
        for(const Particle& p : p_new)
        {
            m_particles.push_back(p);
            m_particleNumber++;
        }
    }
}

void System::applySinks()
{
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
    std::cout << "rest density " << m_fluidParams.restDensity() << std::endl;
    std::cout << "rho avg : " << m_solverParams.averageDensity() << std::endl;
    std::cout << "l : " << l << std::endl;
}

void System::debugFluid()
{
    std::cout << "Particle Number : " << m_particleNumber << std::endl;
    std::cout << "Boundary Number : " << m_boundaryNumber << std::endl;
    std::cout << "Smoothing Radius : " << m_fluidParams.smoothingRadius() << std::endl;
    std::cout << "Speed sound : " << m_fluidParams.soundSpeed() << std::endl;
    std::cout << "Timestep : " << m_solverParams.timeStep() << std::endl;

    std::cout << std::endl;
    m_gridInfo.info();
}

void System::write(const char * filename, std::vector<HReal> data)
{
    std::ofstream outputFile;
    outputFile.open(filename);
    outputFile.precision(16);
    for(unsigned int i=0; i <data.size(); ++i)
    {
        outputFile << data[i] << "\n";
    }
    outputFile.close();
}

void System::write(const char * filename, std::vector<Vec3r > data)
{
    std::ofstream outputFile;
    outputFile.open(filename);
    outputFile.precision(16);
    for(unsigned int i=0; i <data.size(); ++i)
    {
        outputFile << data[i][0] << " " << data[i][1] <<" " <<  data[i][2] << "\n";
    }
    outputFile.close();
}

void System::exportState(const char * baseName)
{
    std::vector< Vec3r > x = getPosition();
    std::vector< Vec3r > v = getVelocity();
    std::vector< HReal > d = getDensity();
    std::vector< HReal > m = getMass();

    std::stringstream ss_padding;
    ss_padding.fill('0');
    ss_padding.width(5);
    ss_padding << m_countExport++;
    std::string padding = ss_padding.str();

    std::stringstream posFilename, velFilename, densFilename, massFilename;
    posFilename << baseName << "/position/position" << padding << ".txt";
    velFilename << baseName << "/velocity/velocity" << padding << ".txt";
    densFilename << baseName << "/density/density" << padding << ".txt";
    massFilename << baseName << "/mass/mass" << padding << ".txt";

    write(  posFilename.str().c_str(), x );
    write(  velFilename.str().c_str(), v );
    write(  densFilename.str().c_str(), d );
    write(  massFilename.str().c_str(), m );
}

HReal System::getTime()
{
    return m_time;
}

HReal & System::getMeanDensityValue()
{
    return m_meanDensity;
}

HReal & System::getDensityFluctuationValue()
{
    return m_densityFluctuation;
}

HReal & System::getRealVolumeValue()
{
    return m_realVolume;
}

const int & System::particleNumber() const
{
    return m_particleNumber;
}

}//namespace hokusai
