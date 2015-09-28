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

    //a_kernel = AkinciKernel();
    p_kernel = MonaghanKernel();
    b_kernel = BoundaryKernel();

    particles = vector<Particle>();
    boundaries = vector<Boundary>();
}

System::System(int resolution)
{
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

    //a_kernel = AkinciKernel();
    p_kernel = MonaghanKernel();
    b_kernel = BoundaryKernel();

    particles = vector<Particle>();
    boundaries = vector<Boundary>();

    setParameters(resolution, 1.0);
}

System::~System(){}

void System::computeRho(int i)
{
    Particle& pi=particles[i];
    pi.rho=0.0;
    for(int j=0; j<(int)particles.size(); ++j)
    {
        Particle& pj=particles[j];
        pi.rho += mass*p_kernel.monaghanValue(pi.x-pj.x);
    }
    for(int j=0; j<(int)boundaries.size(); ++j)
    {
        Boundary& bj = boundaries[j];
        pi.rho += p_kernel.monaghanValue(pi.x-bj.x)*bj.psi;
    }
}

void System::computeNormal(int i)
{
    //Compute normal
    Particle& pi=particles[i];
    Vec2d n(0.0);
    Vec2d gradient(0.0);
    for(int j=0; j< (int)particles.size(); ++j)
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

    for(int j=0; j<(int)particles.size(); ++j)
    {
        computeViscosityForces(i, j);
        computeSurfaceTensionForces(i, j);
    }

    for(int j=0; j<(int)boundaries.size(); ++j)
    {
        computeBoundaryFrictionForces(i, j);
        computeBoundaryAdhesionForces(i, j);
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
    double fdrho=0.0, bdrho=0.0;
    Vec2d gradient = Vec2d::Zero();


    for(int j=0; j<(int)particles.size(); ++j)
    {
        if(i!=j)
        {
            Particle& pj=particles[j];
            p_kernel.monaghanGradient(pi.x-pj.x, gradient);
            Vec2d vij_adv=pi.v_adv-pj.v_adv;
            fdrho+=mass*vij_adv.dot(gradient);
        }
    }

    for(int j=0; j<(int)boundaries.size(); ++j)
    {
        Boundary& bj=boundaries[j];
        Vec2d vb = 0.1*Vec2d::Ones();
        Vec2d v = pi.v_adv - vb; //vb(t+dt)
        p_kernel.monaghanGradient(pi.x-bj.x, gradient);
        bdrho+=bj.psi*v.dot(gradient);
    }

    pi.rho_adv = pi.rho + dt*( fdrho + bdrho );
}

void System::computeSumDijPj(int i)
{
        Particle& pi=particles[i];
        pi.sum_dij = Vec2d::Zero();
        for(int j=0; j<(int)particles.size(); ++j)
        {
            if(i!=j)
            {
                Particle& pj=particles[j];
                Vec2d gradient(0.0);
                p_kernel.monaghanGradient(pi.x-pj.x, gradient);
                pi.sum_dij+=(-mass/pow(pj.rho,2))*pj.p_l*gradient;
            }
        }
        pi.sum_dij *= pow(dt,2);
}

void System::computePressure(int i)
{
        Particle& pi=particles[i];
        double fsum=0.0, bsum=0.0, omega=0.5;

        for(int j=0; j<(int)particles.size(); ++j)
        {
            if(i!=j)
            {
                Particle& pj=particles[j];
                Vec2d gradient_ij(0.0), dji=computeDij(j, i);
                p_kernel.monaghanGradient(pi.x-pj.x, gradient_ij);
                Vec2d aux = pi.sum_dij - (pj.dii_fluid+pj.dii_boundary)*pj.p_l - (pj.sum_dij - dji*pi.p_l);
                fsum+=mass*aux.dot(gradient_ij);
            }
        }

        for(int j=0; j<(int)boundaries.size(); ++j)
        {
            Boundary& bj=boundaries[j];
            Vec2d gradient(0.0), r(0.0); r=pi.x-bj.x;
            p_kernel.monaghanGradient(r, gradient);
            bsum+=bj.psi*pi.sum_dij.dot(gradient);
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
        //Vec2d gradient(0.0);
        pi.f_p = Vec2d::Zero();

        //Fluid Pressure Force
        for(int j=0; j<(int)particles.size(); ++j)
        {
            computeFluidPressureForce(i, j);
        }

        //Boundary Pressure Force [Akinci 2012]
        for(int j=0; j<(int)boundaries.size(); ++j)
        {
            computeBoundaryPressureForce(i, j);
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


void System::computeDii_Boundary(int i)
{
    //    Particle& pi=particles[i];
    //    pi.dii_boundary = Vec2d::Zeros();
    //    for(int& j : pi.boundaryNeighbor)
    //    {
    //        Boundary& bj=boundaries[j];
    //        Vec2d gradient(0.0);
    //        p_kernel.monaghanGradient(pi.x-bj.x, gradient);
    //        pi.dii_boundary+=(-dt*dt*bj.psi/pow(pi.rho,2))*gradient;
    //    }
    //}

    //void System::computeDii_Fluid(int i)
    //{
    //    Particle& pi=particles[i];
    //    pi.dii_fluid = Vec2d::Zeros();
    //    pi.dii_boundary = Vec2d::Zeros();
    //    for(int& j : pi.fluidNeighbor)
    //    {
    //        if(i!=j)
    //        {
    //            Particle& pj=particles[j];
    //            Vec2d gradient(0.0);
    //            p_kernel.monaghanGradient(pi.x-pj.x, gradient);
    //            pi.dii_fluid+=(-dt*dt*mass/pow(pi.rho,2))*gradient;
    //        }
    //    }
}

void System::computeDii(int i)
{
    Particle& pi=particles[i];
    pi.dii_fluid = Vec2d::Zero();
    pi.dii_boundary = Vec2d::Zero();

    for(int j=0; j<(int)particles.size();++j)
    {
        if(i!=j)
        {
            Particle& pj=particles[j];
            Vec2d gradient(0.0);
            p_kernel.monaghanGradient(pi.x-pj.x, gradient);
            pi.dii_fluid+=(-dt*dt*mass/pow(pi.rho,2))*gradient;
        }
    }

    for(int j=0; j<(int)boundaries.size();++j)
    {
        Boundary& bj=boundaries[j];
        Vec2d gradient(0.0);
        p_kernel.monaghanGradient(pi.x-bj.x, gradient);
        pi.dii_boundary+=(-dt*dt*bj.psi/pow(pi.rho,2))*gradient;
    }
}

void System::computeAii( int i)
{
    Particle& pi=particles[i]; pi.aii=0.0;
    for(int j=0; j<(int)particles.size();++j)
    {
        if(i!=j)
        {
            Particle& pj=particles[j];
            Vec2d dji=computeDij(j,i);
            Vec2d gradient_ij(0.0);
            p_kernel.monaghanGradient(pi.x-pj.x, gradient_ij);
            pi.aii+=mass*((pi.dii_fluid+pi.dii_boundary)-dji).dot(gradient_ij);
        }
    }
    for(int j=0; j<(int)boundaries.size();++j)
    {
        Boundary& bj=boundaries[j];
        Vec2d gradient_ij(0.0);
        p_kernel.monaghanGradient(pi.x-bj.x, gradient_ij);
        pi.aii+=bj.psi*(pi.dii_fluid+pi.dii_boundary).dot(gradient_ij);
    }
}

void System::computeBoundaryVolume()
{
    for(int i=0; i<(int)boundaries.size();++i)
    {
        double densityNumber=0.0;
        for(int j=0; j<(int)boundaries.size();++j)
        {
            densityNumber += p_kernel.monaghanValue(boundaries[i].x-boundaries[j].x);
            boundaries[i].psi = restDensity/densityNumber;
        }
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

void System::addBoundaryBox(Vec2d offset, Vec2d scale)
{
    int widthSize = std::floor(scale[0]/h);
    int heightSize = std::floor(scale[1]/h);
    //Bottom+Top
    for(int i=0; i<widthSize; ++i)
    {
        Vec2d bPosition(offset[0]+i*h, offset[1]);
        boundaries.push_back(Boundary(bPosition,Vec2d(0.0),0.0));
        boundaryNumber++;
        Vec2d tPosition(offset[0]+i*h, offset[1]+scale[1]);
        boundaries.push_back(Boundary(tPosition,Vec2d(0.0),0.0));
        boundaryNumber++;

    }
    //Left+Right
    for(int i=0; i<=heightSize; ++i)
    {
        Vec2d lPosition(offset[0],offset[1]+i*h);
        boundaries.push_back(Boundary(lPosition,Vec2d(0.0),0.0));
        boundaryNumber++;
        Vec2d rPosition(offset[0]+scale[0],offset[1]+i*h);
        boundaries.push_back(Boundary(rPosition,Vec2d(0.0),0.0));
        boundaryNumber++;

    }
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

const Vec2d& System::getGravity()
{
    return gravity;
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
    restDensity = 1000;
    mass = (restDensity * volume) / _wishedNumber;
    particlePerCell = 33.8; //better
    h = 0.5*pow( double(volume*particlePerCell) / double(M_PI*_wishedNumber), 1.0/2.0);

    double eta = 0.01;
    double H = 0.1;
    double vf = sqrt( 2*9.81*H );
    cs = vf/(sqrt(eta));

    alpha = 0.1;
    fcohesion = 0.05;
    badhesion = 0.001;
    sigma=1.0;
    boundaryH = h/2.0; //boundaryH must be <= h (neighbor search purpose)

    dt = 1e-10;

    p_kernel = MonaghanKernel( h );
    //a_kernel = AkinciKernel( 2.0*h );
    b_kernel = BoundaryKernel( boundaryH, cs );
}

void System::init()
{
    computeBoundaryVolume();
}

void System::predictAdvection()
{
    for(int i=0; i<particleNumber; ++i)
        computeRho(i);

    for(int i=0; i<particleNumber; ++i)
        computeNormal(i);

    //computeSurfaceParticle();

    for(int i=0; i<particleNumber; ++i)
    {
        computeAdvectionForces(i);
        predictVelocity(i);
        computeDii(i);
    }

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
        for(int i=0; i<particleNumber; ++i)
            computeSumDijPj(i);

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

    for(int i=0; i<particleNumber; ++i)
        computePressureForce(i);

    for(int i=0; i<particleNumber; ++i)
    {
        Particle& pi=particles[i];
        pi.v = pi.v_adv + (dt*pi.f_p)/mass;
        pi.x += dt*pi.v;
    }
}

void System::simulate()
{
    //prepareGrid();
    predictAdvection();
    pressureSolve();
    integration();
    //applySources();
    //applySinks();
    computeStats();
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
