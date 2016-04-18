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

template<typename Solver>
System<Solver>::System() :
    m_particles( ParticleContainerPtr<Particle>() ),
    m_boundaries( BoundaryContainerPtr<Boundary>() )
{
    m_countExport = 0;
    m_countTime = 0;
    m_particleNumber = 0;
    m_boundaryNumber = 0;

    m_volume = 0.0;
    m_restDensity = 0.0;
    m_meanDensity = 0.0;
    m_densityFluctuation = 0.0;
    m_realVolume = 0.0;
    m_mass = 0.0;
    m_h = 0.0;
    m_fcohesion = 0.0;
    m_badhesion = 0.0;
    m_cs = 0.0;
    m_alpha = 0.0;
    m_boundaryH = 0.0;
    m_dt = 0.0;
    m_time = 0.0;
    m_rho_avg_l = 0.0;
    m_maxEta = 1.0;

    m_gravity = Vec3r(0,-9.81,0);

    m_aKernel = AkinciKernel();
    m_pKernel = MonaghanKernel();
    m_bKernel = BoundaryKernel();

    m_gridInfo = GridUtility();
    //m_particles = ParticleContainer<Particle>();
    //m_particles = vector<Particle>();
    //m_boundaries = BoundaryContainer<Boundary>();
    //m_boundaries = std::vector<Boundary>();
}

template<typename Solver>
System<Solver>::System(int resolution) :
    m_particles( ParticleContainerPtr<Particle>() ),
    m_boundaries( BoundaryContainerPtr<Boundary>() )
{
    m_countExport = 0;
    m_countTime = 0;
    m_particleNumber = 0;
    m_boundaryNumber = 0;

    m_volume = 0.0;
    m_restDensity = 0.0;
    m_meanDensity = 0.0;
    m_densityFluctuation = 0.0;
    m_realVolume = 0.0;
    m_mass = 0.0;
    m_h = 0.0;
    m_fcohesion = 0.0;
    m_badhesion = 0.0;
    m_cs = 0.0;
    m_alpha = 0.0;
    m_boundaryH = 0.0;
    m_dt = 0.0;
    m_time = 0.0;
    m_rho_avg_l = 0.0;
    m_maxEta = 1.0;

    m_gravity = Vec3r(0,-9.81,0);

    m_aKernel = AkinciKernel();
    m_pKernel = MonaghanKernel();
    m_bKernel = BoundaryKernel();

    m_gridInfo = GridUtility();
    //m_particles = ParticleContainer<Particle>();
    //m_particles = std::vector<Particle>();
    //m_boundaries = BoundaryContainer<Boundary>();
    //m_boundaries = std::vector<Boundary>();

    setParameters(resolution, 1.0);
}

template<typename Solver>
System<Solver>::~System(){}

template<typename Solver>
void System<Solver>::computeDensity(int i)
{
    Particle& pi=m_particles[i];
    std::vector<int>& fneighbors=pi.fluidNeighbor;
    std::vector<int>& bneighbors=pi.boundaryNeighbor;
    pi.rho=0.0;
    for(int& j : fneighbors)
    {
        Particle& pj=m_particles[j];
        pi.rho += m_mass*m_pKernel.monaghanValue(pi.x-pj.x);
    }
    for(int& j : bneighbors)
    {
        Boundary& bj = m_boundaries[j];
        pi.rho += m_pKernel.monaghanValue(pi.x-bj.x)*bj.psi;
    }
}

template<typename Solver>
void System<Solver>::computeNormal(int i)
{
    //Compute normal
    Particle& pi=m_particles[i];
    std::vector< int > & neighbors = m_particles[i].fluidNeighbor;
    Vec3r n(0.0);
    Vec3r gradient(0.0);
    for(int& j : neighbors)
    {
        if(i!=j)
        {
            Particle& pj = m_particles[j];
            m_pKernel.monaghanGradient(pi.x-pj.x, gradient);
            n += (m_mass/pj.rho)*gradient;
        }
    }
    pi.n = m_h*n;
}

template<typename Solver>
bool System<Solver>::isSurfaceParticle(int i, HReal treshold)
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

template<typename Solver>
std::vector<typename Solver::Particle> System<Solver>::getSurfaceParticle()
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

template<typename Solver>
void System<Solver>::computeAdvectionForces(int i)
{
    Particle& pi=m_particles[i];
    pi.f_adv.fill(0.0);
    for(int& j : pi.fluidNeighbor)
    {
        computeViscosityForces(i, j);
        computeSurfaceTensionForces(i, j);
    }
    for(int& j : pi.boundaryNeighbor)
    {
        computeBoundaryFrictionForces(i, j);
        computeBoundaryAdhesionForces(i, j);
    }
    pi.f_adv+=m_gravity*m_mass;
}

template<typename Solver>
void System<Solver>::predictVelocity(int i)
{
    Particle& pi=m_particles[i];
    pi.v_adv = pi.v + (m_dt/m_mass)*pi.f_adv;
}

template<typename Solver>
void System<Solver>::predictDensity(int i)
{
    Particle& pi=m_particles[i];
    std::vector<int>& fneighbors=pi.fluidNeighbor;
    std::vector<int>& bneighbors=pi.boundaryNeighbor;
    HReal fdrho=0.0, bdrho=0.0;
    Vec3r gradient(0.0);

    for(int& j : fneighbors)
    {
        if(i!=j)
        {
            Particle& pj=m_particles[j];
            m_pKernel.monaghanGradient(pi.x-pj.x, gradient);
            Vec3r vij_adv=pi.v_adv-pj.v_adv;
            fdrho+=m_mass*Vec3r::dotProduct(vij_adv, gradient);
        }
    }

    for(int& j: bneighbors)
    {
        Boundary& bj=m_boundaries[j];
        Vec3r vb(0.1), v(0.0); v = pi.v_adv - vb; //vb(t+dt)
        m_pKernel.monaghanGradient(pi.x-bj.x, gradient);
        bdrho+=bj.psi*Vec3r::dotProduct(v,gradient);
    }

    pi.rho_adv = pi.rho + m_dt*( fdrho + bdrho );
}

template<typename Solver>
void System<Solver>::computeSumDijPj(int i)
{
    Particle& pi=m_particles[i];
    std::vector<int>& fneighbors=pi.fluidNeighbor;
    pi.sum_dij.fill(0.0);
    for(int& j : fneighbors)
    {
        if(i!=j)
        {
            Particle& pj=m_particles[j];
            Vec3r gradient(0.0);
            m_pKernel.monaghanGradient(pi.x-pj.x, gradient);
            pi.sum_dij+=(-m_mass/pow(pj.rho,2))*pj.p_l*gradient;
        }
    }
    pi.sum_dij *= pow(m_dt,2);
}

template<typename Solver>
void System<Solver>::computeViscosityForces(int i, int j)
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

template<typename Solver>
void System<Solver>::computeBoundaryFrictionForces(int i, int j)
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

template<typename Solver>
void System<Solver>::computeSurfaceTensionForces(int i, int j)
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

template<typename Solver>
void System<Solver>::computeBoundaryAdhesionForces(int i, int j)
{
    Particle& pi=m_particles[i];
        Boundary& bj=m_boundaries[j];
        Vec3r xij= pi.x - bj.x;
        HReal l = xij.length();
        pi.f_adv += -(m_badhesion*m_mass*m_boundaries[j].psi*m_aKernel.adhesionValue(l)/l)*xij;
}

template<typename Solver>
Vec3r System<Solver>::computeDij(int i, int j)
{
    Particle& pi=m_particles[i];
    Particle& pj=m_particles[j];
    Vec3r gradient(0.0);
    m_pKernel.monaghanGradient(pi.x-pj.x, gradient);
    Vec3r d=-(m_dt*m_dt*m_mass)/pow(pj.rho,2)*gradient;
    return d;
}

template<typename Solver>
void System<Solver>::computePressure(int i)
{
    Particle& pi=m_particles[i];
    std::vector<int>& fneighbors=pi.fluidNeighbor, bneighbors=pi.boundaryNeighbor;
    HReal fsum=0.0, bsum=0.0, omega=0.5;

    for(int& j : fneighbors)
    {
        if(i!=j)
        {
            Particle& pj=m_particles[j];
            Vec3r gradient_ij(0.0), dji=computeDij(j, i);
            m_pKernel.monaghanGradient(pi.x-pj.x, gradient_ij);
            Vec3r aux = pi.sum_dij - (pj.dii_fluid+pj.dii_boundary)*pj.p_l - (pj.sum_dij - dji*pi.p_l);
            fsum+=m_mass*Vec3r::dotProduct(aux, gradient_ij);
        }
    }

    for(int& j : bneighbors)
    {
        Boundary& bj=m_boundaries[j];
        Vec3r gradient(0.0), r(0.0); r=pi.x-bj.x;
        m_pKernel.monaghanGradient(r, gradient);
        bsum+=bj.psi*Vec3r::dotProduct(pi.sum_dij,gradient);
    }

    HReal previousPl = pi.p_l;
    pi.rho_corr = pi.rho_adv + fsum + bsum;
    if(std::abs(pi.aii)>std::numeric_limits<HReal>::epsilon())
        pi.p_l = (1-omega)*previousPl + (omega/pi.aii)*(m_restDensity - pi.rho_corr);
    else
        pi.p_l = 0.0;
    pi.p = std::max(pi.p_l,0.0);
    pi.p_l = pi.p;
    pi.rho_corr += pi.aii*previousPl;
}

template<typename Solver>
void System<Solver>::computePressureForce(int i)
{
    Particle& pi=m_particles[i];
    pi.f_p.fill(0.0);
    std::vector<int>& fneighbors=pi.fluidNeighbor;
    std::vector<int>& bneighbors=pi.boundaryNeighbor;

    //Fluid Pressure Force
    for(int& j : fneighbors)
    {
        computeFluidPressureForce(i, j);
    }

    //Boundary Pressure Force [Akinci 2012]
    for(int& j : bneighbors )
    {
        computeBoundaryPressureForce(i, j);
    }
}

template<typename Solver>
void System<Solver>::computeFluidPressureForce(int i, int j)
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

template<typename Solver>
void System<Solver>::computeBoundaryPressureForce(int i, int j)
{
    Vec3r gradient(0.0);
    Particle& pi=m_particles[i];
    Boundary& bj=m_boundaries[j];
        m_pKernel.monaghanGradient(pi.x-bj.x, gradient);
        pi.f_p += -m_mass*bj.psi*( pi.p/pow(pi.rho,2) ) * gradient;
}

template<typename Solver>
void System<Solver>::initializePressure(int i)
{
    Particle& pi=m_particles[i];
    pi.p_l=0.5*pi.p;
}

template<typename Solver>
void System<Solver>::computeError()
{
    m_rho_avg_l=0.0;
    for(int i=0; i<m_particleNumber; ++i)
        m_rho_avg_l+=m_particles[i].rho_corr;
    m_rho_avg_l /= m_particleNumber;
}


template<typename Solver>
void System<Solver>::computeDii_Boundary(int i)
{
    Particle& pi=m_particles[i];
    pi.dii_boundary.fill(0.0);
    for(int& j : pi.boundaryNeighbor)
    {
        Boundary& bj=m_boundaries[j];
        Vec3r gradient(0.0);
        m_pKernel.monaghanGradient(pi.x-bj.x, gradient);
        pi.dii_boundary+=(-m_dt*m_dt*bj.psi/pow(pi.rho,2))*gradient;
    }
}

template<typename Solver>
void System<Solver>::computeDii_Fluid(int i)
{
    Particle& pi=m_particles[i];
    pi.dii_fluid.fill(0.0);
    pi.dii_boundary.fill(0.0);
    for(int& j : pi.fluidNeighbor)
    {
        if(i!=j)
        {
            Particle& pj=m_particles[j];
            Vec3r gradient(0.0);
            m_pKernel.monaghanGradient(pi.x-pj.x, gradient);
            pi.dii_fluid+=(-m_dt*m_dt*m_mass/pow(pi.rho,2))*gradient;
        }
    }
}

template<typename Solver>
void System<Solver>::computeDii(int i)
{
    Particle& pi=m_particles[i];
    pi.dii_fluid.fill(0.0);
    pi.dii_boundary.fill(0.0);
    for(int& j : pi.fluidNeighbor)
    {
        if(i!=j)
        {
            Particle& pj=m_particles[j];
            Vec3r gradient(0.0);
            m_pKernel.monaghanGradient(pi.x-pj.x, gradient);
            pi.dii_fluid+=(-m_dt*m_dt*m_mass/pow(pi.rho,2))*gradient;
        }
    }
    for(int& j : pi.boundaryNeighbor)
    {
        Boundary& bj=m_boundaries[j];
        Vec3r gradient(0.0);
        m_pKernel.monaghanGradient(pi.x-bj.x, gradient);
        pi.dii_boundary+=(-m_dt*m_dt*bj.psi/pow(pi.rho,2))*gradient;
    }
}

template<typename Solver>
void System<Solver>::computeAii( int i)
{
    Particle& pi=m_particles[i]; pi.aii=0.0;
    for(int& j : pi.fluidNeighbor)
    {
        if(i!=j)
        {
            Particle& pj=m_particles[j];
            Vec3r dji=computeDij(j,i);
            Vec3r gradient_ij(0.0);
            m_pKernel.monaghanGradient(pi.x-pj.x, gradient_ij);
            pi.aii+=m_mass*Vec3r::dotProduct((pi.dii_fluid+pi.dii_boundary)-dji,gradient_ij);
        }
    }
    for(int& j : pi.boundaryNeighbor)
    {
        Boundary& bj=m_boundaries[j];
        Vec3r gradient_ij(0.0);
        m_pKernel.monaghanGradient(pi.x-bj.x, gradient_ij);
        pi.aii+=bj.psi*Vec3r::dotProduct(pi.dii_fluid+pi.dii_boundary,gradient_ij);
    }
}

template<typename Solver>
void System<Solver>::getNearestNeighbor(const int i, const HReal radius)
{
    Particle& p = m_particles[i];
    p.fluidNeighbor.clear();
    p.boundaryNeighbor.clear();

    std::vector<int> neighborCell;
    m_gridInfo.get27Neighbors(neighborCell, p.x, radius);

    for(size_t i=0; i<neighborCell.size(); ++i)
    {
        std::vector<int>& bNeighborCell = m_boundaryGrid[neighborCell[i]];
        for(size_t j=0; j<bNeighborCell.size(); ++j)
        {
            int bParticleId = m_boundaryGrid[neighborCell[i]][j];
            Boundary& bParticle = m_boundaries[bParticleId];
            Vec3r d = bParticle.x-p.x;
            if( d.lengthSquared()<radius*radius )
                p.boundaryNeighbor.push_back(bParticleId);
        }

        std::vector<int>& fNeighborCell = m_fluidGrid[neighborCell[i]];
        for(size_t j=0; j< fNeighborCell.size(); ++j)
        {
            int fParticleId = m_fluidGrid[neighborCell[i]][j];
            Particle& fParticle = m_particles[fParticleId];
            Vec3r d = fParticle.x-p.x;
            if( d.lengthSquared()<radius*radius)
                p.fluidNeighbor.push_back(fParticleId);
        }
    }
}

template<typename Solver>
ParticleContainer<typename Solver::Particle>& System<Solver>::getParticles()
{
    return *(m_particles.get());
}

template<typename Solver>
void System<Solver>::getNearestNeighbor(std::vector< int >& neighbor, const std::vector< std::vector<int> >& grid, const Vec3r &x)
{
    std::vector<int> neighborCell;
    m_gridInfo.get27Neighbors(neighborCell, x, m_gridInfo.spacing());
    neighbor.clear();
    for(size_t i=0; i<neighborCell.size(); ++i)
    {
        for(size_t j=0; j<grid[neighborCell[i]].size(); ++j)
        {
            neighbor.push_back(grid[neighborCell[i]][j]);
        }
    }
}

template<typename Solver>
void System<Solver>::computeBoundaryVolume()
{
    for(int i=0; i<m_boundaryNumber; ++i)
    {
        HReal densityNumber=0.0;
        std::vector<int> neighbors;
        getNearestNeighbor(neighbors, m_boundaryGrid, m_boundaries[i].x);
        for(int& j : neighbors)
            densityNumber += m_pKernel.monaghanValue(m_boundaries[i].x-m_boundaries[j].x);
        m_boundaries[i].psi = m_restDensity/densityNumber;
    }
}

template<typename Solver>
void System<Solver>::computeMeanDensity()
{
    m_meanDensity=0.0;
    for(int i=0; i<m_particleNumber; ++i)
    {
        m_meanDensity+=m_particles[i].rho;
    }
    m_meanDensity/=m_particleNumber;
}

template<typename Solver>
void System<Solver>::computeDensityFluctuation()
{
    m_densityFluctuation=m_meanDensity-m_restDensity;
}

template<typename Solver>
void System<Solver>::computeVolume()
{
    m_realVolume=0.0;
    for(int i=0; i<m_particleNumber; ++i)
    {
        m_realVolume += m_mass/m_particles[i].rho;
    }
}

template<typename Solver>
void System<Solver>::setGravity(const Vec3r& _gravity)
{
    m_gravity = _gravity;
}

template<typename Solver>
const Vec3r& System<Solver>::getGravity()
{
    return m_gravity;
}

template<typename Solver>
void System<Solver>::setParameters( int _wishedNumber, HReal _volume )
{
    m_time = 0.0;
    m_countTime = 0.0;
    m_countExport = 0;
    m_particleNumber = 0;
    m_boundaryNumber = 0;
    m_volume = _volume;

    m_maxEta=1.0;
    m_restDensity = 1000;
    m_mass = (m_restDensity * m_volume) / _wishedNumber;
    m_particlePerCell = 33.8; //better
    m_h = 0.5*pow( HReal(3*m_volume*m_particlePerCell) / HReal(4*M_PI*_wishedNumber), 1.0/3.0);

    HReal eta = 0.01;
    HReal H = 0.1;
    HReal vf = sqrt( 2*9.81*H );
    m_cs = vf/(sqrt(eta));

    m_alpha = 0.1;
    m_fcohesion = 0.05;
    m_badhesion = 0.001;
    m_sigma=1.0;
    m_boundaryH = m_h/2.0; //boundaryH must be <= h (neighbor search purpose)

    m_dt = 0.004;

    m_pKernel = MonaghanKernel( m_h );
    m_aKernel = AkinciKernel( 2.0*m_h );
    m_bKernel = BoundaryKernel( m_boundaryH, m_cs );
}

template<typename Solver>
void System<Solver>::init()
{
    mortonSort();

    m_boundaryGrid.resize(m_gridInfo.size());
    for(size_t i=0; i<m_boundaries.size(); ++i)
    {
        int id = m_gridInfo.cellId(m_boundaries[i].x);
        if(m_gridInfo.isInside(id))
            m_boundaryGrid[id].push_back(i);
    }

    m_fluidGrid.resize(m_gridInfo.size());

    //Init simulation values
    computeBoundaryVolume();
    prepareGrid();

    for(size_t i=0; i<m_particles.size(); ++i)
    {
        m_particles[i].isSurface = true;
    }

    debugFluid();
}

template<typename Solver>
void System<Solver>::addBoundaryBox(const Vec3r& offset, const Vec3r& scale)
{
    int epsilon = 0;
    int widthSize = floor(scale[0]/m_h);
    int heightSize = floor(scale[1]/m_h);
    int depthSize = floor(scale[2]/m_h);

    //ZX plane - bottom
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec3r position(i*m_h, offset[1], j*m_h);
            m_boundaries.add(Boundary(position,Vec3r(0.0),0.0));
            m_boundaryNumber++;
        }
    }

    //ZX plane - top
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec3r position(i*m_h, offset[1]+scale[1], j*m_h);
            m_boundaries.add(Boundary(position,Vec3r(0.0),0.0));
            m_boundaryNumber++;
        }
    }

    //XY plane - back
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= heightSize+epsilon; ++j)
        {
            Vec3r position(i*m_h, j*m_h, offset[2]);
            m_boundaries.add(Boundary(position,Vec3r(0.0),0.0));
            m_boundaryNumber++;
        }
    }

    //XY plane - front
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= heightSize-epsilon; ++j)
        {
            Vec3r position(i*m_h, j*m_h, offset[2]+scale[2]);
            m_boundaries.add(Boundary(position,Vec3r(0.0),0.0));
            m_boundaryNumber++;
        }
    }

    //YZ plane - left
    for(int i = -epsilon; i <= heightSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec3r position(offset[0], i*m_h, j*m_h);
            m_boundaries.add(Boundary(position,Vec3r(0.0),0.0));
            m_boundaryNumber++;
        }
    }

    //YZ plane - right
    for(int i = -epsilon; i <= heightSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec3r position(offset[0]+scale[0], i*m_h, j*m_h);
            m_boundaries.add(Boundary(position,Vec3r(0.0),0.0));
            m_boundaryNumber++;
        }
    }

    m_gridInfo.update(offset-Vec3r(2.0*m_h), scale+Vec3r(4.0*m_h), 2.0*m_h);
}

template<typename Solver>
void System<Solver>::addBoundarySphere(const Vec3r& offset, const HReal& radius)
{
    std::vector<Vec3r> samples = getSphereSampling(offset, radius, m_h, m_h);
    for(size_t i=0; i<samples.size(); ++i)
    {
        m_boundaries.add(Boundary(samples[i],Vec3r(0.0),0.0));
        m_boundaryNumber++;
    }
}

template<typename Solver>
void System<Solver>::addBoundaryHemiSphere(const Vec3r& offset, const HReal& radius)
{
    std::vector<Vec3r> samples = getHemiSphereSampling(offset, radius, m_h, m_h);
    for(size_t i=0; i<samples.size(); ++i)
    {
        m_boundaries.add(Boundary(samples[i],Vec3r(0.0),0.0));
        m_boundaryNumber++;
    }
}

template<typename Solver>
void System<Solver>::addBoundaryDisk(const Vec3r& offset, const HReal& radius)
{
    std::vector<Vec3r> samples = getDiskSampling(offset, radius, m_h);
    for(size_t i=0; i<samples.size(); ++i)
    {
        m_boundaries.add(Boundary(samples[i],Vec3r(0.0),0.0));
        m_boundaryNumber++;
    }
}


template<typename Solver>
void System<Solver>::translateBoundaries(const Vec3r& t)
{
    for(size_t i=0; i<m_boundaries.size(); ++i)
    {
        m_boundaries[i].x += t;
    }
}

template<typename Solver>
void System<Solver>::translateParticles(const Vec3r& t)
{
    for(size_t i=0; i<m_particles.size(); ++i)
    {
        m_particles[i].x += t;
    }
}

template<typename Solver>
void System<Solver>::addParticleSphere(const Vec3r& centre, const HReal radius)
{
    Vec3r scale(2.0*radius, 2.0*radius, 2.0*radius);
    Vec3r offset = centre - Vec3r(radius, radius, radius);
    GridUtility grid(offset, scale, m_h);

    for(int i=0; i<grid.size(); ++i)
    {
        Vec3r _x = grid.gridToWorld(i);
        _x+=grid.spacing()/2.0;
        Vec3r _v(0,0,0);
        HReal l2 = (centre-_x).lengthSquared();
        if(l2<=(radius*radius))
        {
            m_particles.add( Particle(_x,_v) );
            m_particleNumber++;
        }
    }
}

template<typename Solver>
void System<Solver>::addParticleSource(const ParticleSource<Particle> &s)
{
    m_pSources.push_back(s);
}

template<typename Solver>
void System<Solver>::addParticleBox(const Vec3r& offset, const Vec3r& scale)
{
    int widthSize = floor(scale[0]/m_h);
    int heightSize = floor(scale[1]/m_h);
    int depthSize = floor(scale[2]/m_h);

    for(int i=0; i<widthSize; ++i)
    {
        for(int j=0; j<heightSize; ++j)
        {
            for(int k=0; k<depthSize; ++k)
            {
                Vec3r _x = offset + Vec3r(i*m_h,j*m_h,k*m_h);
                Vec3r _v(0,0,0);
                m_particles.add( Particle(_x,_v) );
                m_particleNumber++;
            }
        }
    }
}

template<typename Solver>
void System<Solver>::addParticleBox(HReal width, HReal height, HReal depth, HReal spacing)
{
    int widthSize = floor(width/spacing);
    int heightSize = floor(height/spacing);
    int depthSize = floor(depth/spacing);

    for(int i=0; i<widthSize; ++i)
    {
        for(int j=0; j<heightSize; ++j)
        {
            for(int k=0; k<depthSize; ++k)
            {
                Vec3r _x(i*spacing,j*spacing,k*spacing);
                Vec3r _v(0,0,0);
                m_particles.add( Particle(_x,_v) );
                m_particleNumber++;
            }
        }
    }
}


template<typename Solver>
void System<Solver>::createParticleVolume(Vec3r& pos, HReal width, HReal /*height*/, HReal depth, HReal spacing, int particleMax)
{
    int widthSize = floor( width/spacing );
    int depthSize = floor( depth/spacing );
    int count = 0;
    int j = 0;
    while(count < particleMax)
    {
        for(int i = 0; i <= widthSize; ++i)
        {
            for(int k = 0; k <= depthSize; ++k)
            {
                if(count < particleMax)
                {
                    Vec3r _x(pos[0]+i*spacing,pos[1]+j*spacing,pos[2]+k*spacing);
                    Vec3r _v(0,0,0);
                    m_particles.add( Particle(_x,_v) );
                    m_particleNumber++;
                }
                count++;
            }
        }
        j++;
    }
}

template<typename Solver>
void System<Solver>::addFluidParticle(const Vec3r& x, const Vec3r& v)
{
    m_particles.add( Particle(x,v) );
    m_particleNumber++;
}

template<typename Solver>
void System<Solver>::addBoundaryParticle(const Vec3r& x, const Vec3r& v)
{
    m_boundaries.add( Boundary(x,v) );
    m_boundaryNumber++;
}

template<typename Solver>
void System<Solver>::addBoundaryMesh(const char* filename)
{
    TriMesh mesh(filename);
    std::vector<Vec3r> samples;
    AkinciMeshSampling(mesh, m_h/2.0, samples);
    Vec3r minBB(std::numeric_limits<HReal>::max()), maxBB(-std::numeric_limits<HReal>::max());
    for(size_t i=0; i<samples.size(); ++i)
    {
        for(int j=0; j<3; ++j)
        {
            minBB[j] = std::min(samples[i][j], minBB[j]);
            maxBB[j] = std::max(samples[i][j], maxBB[j]);
        }
        m_boundaries.add(Boundary(samples[i],Vec3r(0.0),0.0));
        m_boundaryNumber++;
    }
    Vec3r offset = minBB;
    Vec3r scale = maxBB-minBB;
    m_gridInfo.update(offset-Vec3r(2.0*m_h), scale+Vec3r(4.0*m_h), 2.0*m_h);
}

bool pairCompare( const std::pair<int,int>& e1, const std::pair<int,int>& e2 )
{
    return (e1.second < e2.second);
}

template<typename Solver>
void System<Solver>::mortonSort()
{
    std::vector< std::pair<int, int> > particleZindex;

    //Fill particleZindex with particle index and Z-index
    for(int i = 0; i < m_particleNumber; ++i)
    {
        Vec3i _gridIndex = m_gridInfo.worldToGrid(m_particles[i].x);
        std::array<int,3> gridIndex;
        for(int j=0; j<3; ++j)
            gridIndex[j] = _gridIndex[j];
        int zindex = mortonNumber( gridIndex );
        std::pair<int,int> paire(i, zindex);
        particleZindex.push_back(paire);
    }

    //Sort according to z-index
    std::sort( particleZindex.begin(), particleZindex.end(), pairCompare );

    //Move particles according to z-index
    ParticleContainerPtr<Particle> oldParticles(m_particles);

    //ParticleContainer<Particle> oldParticles = m_particles;
    //vector< Particle > oldParticles = m_particles;

    //HReal min = particleZindex[0].second;
    //HReal max = particleZindex[particleNumber-1].second;

    for(int i = 0; i < m_particleNumber; ++i)
    {
        std::pair<int,int>& paire = particleZindex[i];
        if(i != paire.first)
        {
            m_particles[i] = oldParticles[paire.first];
        }
        //HReal color = (paire.second-min)/(max-min);
        //particles[i].c[0] = color;
        //particles[i].c[1] = 0;
        //particles[i].c[2] = 0;
    }
}

template<typename Solver>
void System<Solver>::computeSurfaceParticle()
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

template<typename Solver>
void System<Solver>::prepareGrid()
{
    if( m_countTime%100 == 0 )
        mortonSort();

    for(size_t i=0; i<m_fluidGrid.size(); ++i)
        m_fluidGrid[i].clear();

    for(size_t i=0; i<m_particles.size(); ++i)
    {
        int id = m_gridInfo.cellId(m_particles[i].x);
        if(m_gridInfo.isInside(id))
            m_fluidGrid[id].push_back(i);
    }

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < m_particleNumber; ++i)
        getNearestNeighbor(i, 2.0*m_h);
}

template<typename Solver>
void System<Solver>::predictAdvection()
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
        computeAdvectionForces(i);
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

template<typename Solver>
void System<Solver>::pressureSolve()
{
    int l=0; m_rho_avg_l = 0.0;

    while( ( (m_rho_avg_l-m_restDensity)>m_maxEta ) || (l<2) )
    {
#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
        for(int i=0; i<m_particleNumber; ++i)
            computeSumDijPj(i);

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

template<typename Solver>
void System<Solver>::integration()
{
    m_countTime++; m_time+=m_dt;

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<m_particleNumber; ++i)
        computePressureForce(i);

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<m_particleNumber; ++i)
    {
        Particle& pi=m_particles[i];
        pi.v = pi.v_adv + (m_dt*pi.f_p)/m_mass;
        pi.x += m_dt*pi.v;
    }
}

template<typename Solver>
void System<Solver>::computeSimulationStep()
{
    prepareGrid();
    predictAdvection();
    pressureSolve();
    integration();

    applySources();
    applySinks();
    computeStats();
}

template<typename Solver>
void System<Solver>::applySources()
{
    for(ParticleSource<Particle>& s : m_pSources)
    {
        std::vector<Particle> p_new = s.apply(this->m_time);
        for(const Particle& p : p_new)
        {
            m_particles.add(p);
            m_particleNumber++;
        }
    }
}

template<typename Solver>
void System<Solver>::applySinks()
{
}

template<typename Solver>
void System<Solver>::computeStats()
{
    computeMeanDensity();
    computeVolume();
    computeDensityFluctuation();
}

template<typename Solver>
void System<Solver>::debugIteration(int l)
{
    std::cout.precision(10);
    std::cout << "rest density " << m_restDensity << std::endl;
    std::cout << "rho avg : " << m_rho_avg_l << std::endl;
    std::cout << "l : " << l << std::endl;
}

template<typename Solver>
void System<Solver>::debugFluid()
{
    std::cout << "Particle Number : " << m_particleNumber << std::endl;
    std::cout << "Boundary Number : " << m_boundaryNumber << std::endl;
    std::cout << "Smoothing Radius : " << m_h << std::endl;
    std::cout << "Radius : " << pow( 3.0*m_volume/(4*M_PI*m_particleNumber), 1.0/3.0 ) << std::endl;
    std::cout << "Speed sound : " << m_cs << std::endl;
    std::cout << "Timestep : " << m_dt << std::endl;

    std::cout << std::endl;
    m_gridInfo.info();
}

template<typename Solver>
void System<Solver>::write(const char * filename, std::vector<HReal> data)
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

template<typename Solver>
void System<Solver>::write(const char * filename, std::vector<Vec3r > data)
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

template<typename Solver>
void System<Solver>::exportState(const char * baseName)
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

}//namespace hokusai
