#include "../include/hokusai/system.hpp"
#include "../include/hokusai/rasterizer.hpp"
#include "../include/hokusai/gridUtility.hpp"
#include "../include/hokusai/HSL2RGB.hpp"
#include <algorithm>
#include <assert.h>
#include <omp.h>
#include <fstream>
#include <sstream>

//-----------------------------------------------------
//         System : creation and simulation functions
//-----------------------------------------------------

System::System()
{
    countExport = 0;
    countTime = 0;
    particleNumber = 0;
    boundaryNumber = 0;

    gridChange = 0.0;
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
    scene_radius = 0.0;
    rho_avg_l = 0.0;
    maxEta = 1.0;

    scene_center = Vec(0.0);
    gravity = Vec(0,-9.81,0);
    minBoundary = Vec(0,0,0);
    maxBoundary = Vec(0,0,0);

    a_kernel = AkinciKernel();
    p_kernel = MonaghanKernel();
    b_kernel = BoundaryKernel();

    particles = vector<Particle>();
    boundaries = vector<Boundary>();
    grid = Grid3D();
    visuBox = Box();
}

System::System(int resolution)
{
    countExport = 0;
    countTime = 0;
    particleNumber = 0;
    boundaryNumber = 0;

    gridChange = 0.0;
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
    scene_radius = 0.0;
    rho_avg_l = 0.0;
    maxEta = 1.0;

    scene_center = Vec(0.0);
    gravity = Vec(0,-9.81,0);
    minBoundary = Vec(0,0,0);
    maxBoundary = Vec(0,0,0);

    a_kernel = AkinciKernel();
    p_kernel = MonaghanKernel();
    b_kernel = BoundaryKernel();

    particles = vector<Particle>();
    boundaries = vector<Boundary>();
    grid = Grid3D();
    visuBox = Box();

    setParameters(resolution, 1.0);
}

System::~System(){}

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
    Vec n(0.0);
    Vec gradient(0.0);
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

bool System::isSurfaceParticle(int i, double treshold)
{
    double n_length = particles[i].n.lengthSquared();
    if( n_length > treshold )
        return true;
    else
        return false;
}

vector<Particle> System::getSurfaceParticle()
{
#pragma omp parallel for
    for(int i=0; i<particleNumber; ++i)
        computeRho(i);

#pragma omp parallel for
    for(int i=0; i<particleNumber; ++i)
        computeNormal(i);

    double treshold = 0.05;
    vector<Particle> surfaceParticles;
    for(int i=0; i<particleNumber; ++i)
    {
        if( isSurfaceParticle(i, treshold) )
            surfaceParticles.push_back(particles[i]);
    }
    return surfaceParticles;
}

void System::computeAdvectionForces(int i)
{
    Particle& pi=particles[i];
    pi.f_adv.setAllValue(0.0);
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
    //pi.f_adv[1]+=-9.81*mass;
    pi.f_adv+=gravity*mass;
}

void System::computeBoundaryFrictionForces(int i, int j)
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
        double sigma=1.0;
        double nu = (sigma*h*cs)/(2.0*pi.rho);
        double Pij = -nu * ( std::min(dotVijRij,0.0) / (xij.lengthSquared() + epsilon*h*h) );
        p_kernel.monaghanGradient(xij, gradient);
        pi.f_adv += -mass*bj.psi*Pij*gradient;
    }
}

void System::computeBoundaryAdhesionForces(int i, int j)
{
    Particle& pi=particles[i];
    Boundary& bj=boundaries[j];
    Vec xij= pi.x - bj.x;
    double l = xij.length();
    pi.f_adv += -(badhesion*mass*boundaries[j].psi*a_kernel.adhesionValue(l)/l)*xij;
}

void System::computeViscosityForces(int i, int j)
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

void System::computeSurfaceTensionForces(int i, int j)
{
    if(i!=j)
    {
        Particle& pi=particles[i];
        Particle& pj=particles[j];
        Vec r = pi.x - pj.x;
        double kij = 2.0*restDensity/(pi.rho+pj.rho);
        double l = r.length();
        Vec cohesionForce = -(fcohesion*mass*mass*a_kernel.cohesionValue(l)/l) * r;
        Vec nij = pi.n-pj.n;
        Vec curvatureForce = -fcohesion*mass*nij;
        pi.f_adv += kij*(cohesionForce+curvatureForce);
    }
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
    Vec gradient(0.0);

    for(int& j : fneighbors)
    {
        if(i!=j)
        {
            Particle& pj=particles[j];
            p_kernel.monaghanGradient(pi.x-pj.x, gradient);
            Vec vij_adv=pi.v_adv-pj.v_adv;
            fdrho+=mass*Vec::dotProduct(vij_adv, gradient);
        }
    }

    for(int& j: bneighbors)
    {
        Boundary& bj=boundaries[j];
        Vec vb(0.1), v(0.0); v = pi.v_adv - vb; //vb(t+dt)
        p_kernel.monaghanGradient(pi.x-bj.x, gradient);
        bdrho+=bj.psi*Vec::dotProduct(v,gradient);
    }

    pi.rho_adv = pi.rho + dt*( fdrho + bdrho );
}

Vec System::computeDij(int i, int j)
{
    Particle& pi=particles[i];
    Particle& pj=particles[j];
    Vec gradient(0.0);
    p_kernel.monaghanGradient(pi.x-pj.x, gradient);
    Vec d=-(dt*dt*mass)/pow(pj.rho,2)*gradient;
    return d;
}

void System::computeSumDijPj(int i)
{
    Particle& pi=particles[i];
    vector<int>& fneighbors=pi.fluidNeighbor;
    pi.sum_dij.setAllValue(0.0);
    for(int& j : fneighbors)
    {
        if(i!=j)
        {
            Particle& pj=particles[j];
            Vec gradient(0.0);
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

    for(int& j : fneighbors)
    {
        if(i!=j)
        {
            Particle& pj=particles[j];
            Vec gradient_ij(0.0), dji=computeDij(j, i);
            p_kernel.monaghanGradient(pi.x-pj.x, gradient_ij);
            Vec aux = pi.sum_dij - (pj.dii_fluid+pj.dii_boundary)*pj.p_l - (pj.sum_dij - dji*pi.p_l);
            fsum+=mass*Vec::dotProduct(aux, gradient_ij);
        }
    }

    for(int& j : bneighbors)
    {
        Boundary& bj=boundaries[j];
        Vec gradient(0.0), r(0.0); r=pi.x-bj.x;
        p_kernel.monaghanGradient(r, gradient);
        bsum+=bj.psi*Vec::dotProduct(pi.sum_dij,gradient);
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
    Vec gradient(0.0);
    pi.f_p.setAllValue(0.0);
    vector<int>& fneighbors=pi.fluidNeighbor;
    vector<int>& bneighbors=pi.boundaryNeighbor;

    //Fluid Pressure Force
    for(int& j : fneighbors)
    {
        Particle& pj=particles[j];
        p_kernel.monaghanGradient(pi.x-pj.x, gradient);
        if( i!=j )
        {
            pi.f_p += -mass*mass*( pi.p/pow(pi.rho,2) + pj.p/pow(pj.rho,2) ) * gradient;
        }
    }

    //Boundary Pressure Force [Akinci 2012]
    for(int& j : bneighbors )
    {
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


void System::computeDii_Boundary(int i)
{
    Particle& pi=particles[i];
    pi.dii_boundary.setAllValue(0.0);
    for(int& j : pi.boundaryNeighbor)
    {
        Boundary& bj=boundaries[j];
        Vec gradient(0.0);
        p_kernel.monaghanGradient(pi.x-bj.x, gradient);
        pi.dii_boundary+=(-dt*dt*bj.psi/pow(pi.rho,2))*gradient;
    }
}

void System::computeDii_Fluid(int i)
{
    Particle& pi=particles[i];
    pi.dii_fluid.setAllValue(0.0);
    pi.dii_boundary.setAllValue(0.0);
    for(int& j : pi.fluidNeighbor)
    {
        if(i!=j)
        {
            Particle& pj=particles[j];
            Vec gradient(0.0);
            p_kernel.monaghanGradient(pi.x-pj.x, gradient);
            pi.dii_fluid+=(-dt*dt*mass/pow(pi.rho,2))*gradient;
        }
    }
}

void System::computeDii(int i)
{
    Particle& pi=particles[i];
    pi.dii_fluid.setAllValue(0.0);
    pi.dii_boundary.setAllValue(0.0);
    for(int& j : pi.fluidNeighbor)
    {
        if(i!=j)
        {
            Particle& pj=particles[j];
            Vec gradient(0.0);
            p_kernel.monaghanGradient(pi.x-pj.x, gradient);
            pi.dii_fluid+=(-dt*dt*mass/pow(pi.rho,2))*gradient;
        }
    }
    for(int& j : pi.boundaryNeighbor)
    {
        Boundary& bj=boundaries[j];
        Vec gradient(0.0);
        p_kernel.monaghanGradient(pi.x-bj.x, gradient);
        pi.dii_boundary+=(-dt*dt*bj.psi/pow(pi.rho,2))*gradient;
    }
}

void System::computeAii( int i)
{
    Particle& pi=particles[i]; pi.aii=0.0;
    for(int& j : pi.fluidNeighbor)
    {
        if(i!=j)
        {
            Particle& pj=particles[j];
            Vec dji=computeDij(j,i);
            Vec gradient_ij(0.0);
            p_kernel.monaghanGradient(pi.x-pj.x, gradient_ij);
            pi.aii+=mass*Vec::dotProduct((pi.dii_fluid+pi.dii_boundary)-dji,gradient_ij);
        }
    }
    for(int& j : pi.boundaryNeighbor)
    {
        Boundary& bj=boundaries[j];
        Vec gradient_ij(0.0);
        p_kernel.monaghanGradient(pi.x-bj.x, gradient_ij);
        pi.aii+=bj.psi*Vec::dotProduct(pi.dii_fluid+pi.dii_boundary,gradient_ij);
    }
}

void System::computeBoundaryVolume()
{
    for(int i=0; i<boundaryNumber; ++i)
    {
        double densityNumber=0.0;
        vector<int> neighbors = grid.getNearestBoundaryNeighbor(boundaries[i].x);
        for(int& j : neighbors)
            densityNumber += p_kernel.monaghanValue(boundaries[i].x-boundaries[j].x);
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

void System::setGravity(const Vec& _gravity)
{
    gravity = _gravity;
}

const Vec& System::getGravity()
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
    float particleAverage = 33.8; //better
    h = 0.5*pow( double(3*volume*particleAverage) / double(4*M_PI*_wishedNumber), 1.0/3.0);

    double eta = 0.01;
    double H = 0.1;
    double vf = sqrt( 2*9.81*H );
    cs = vf/(sqrt(eta));

    alpha = 0.1;
    fcohesion = 0.05;
    badhesion = 0.001;
    grid.cellSize = 2.0*h;
    boundaryH = h/2.0; //boundaryH must be <= h (neighbor search purpose)

    dt = 0.004;

    p_kernel = MonaghanKernel( h );
    a_kernel = AkinciKernel( 2.0*h );
    b_kernel = BoundaryKernel( boundaryH, cs );
}

void System::init()
{
    //Init grid
    grid.minCorner = minBoundary - 20*h;
    grid.maxCorner = maxBoundary + 20*h;
    grid.init();
    for( int i = 0; i < boundaryNumber; ++i )
        grid.addBoundary(i, boundaries[i].x);
    for(int i = 0; i < particleNumber; ++i)
        grid.addFluid(i, particles[i]);
    mortonSort();

    //Init simulation values
    computeBoundaryVolume();
    prepareGrid(); 

    for(int i=0; i<particleNumber; ++i)
    {
        computeRho(i);
        particles[i].rho = restDensity;
    }

    for(int i=0; i<particleNumber; ++i)
    {
        computeDii(i);
    }

    for(int i=0; i<particleNumber; ++i)
        computeAii(i);

    computeSceneParam();
    debugFluid();
}

void System::createDamBreak(const Vec& offset, const Vec& dimension, const int pNumber)
{
    double volume = dimension[0]*dimension[1]*dimension[2];
    setParameters(pNumber, volume);
    Vec minOffset = offset;
    addParticleBox(minOffset, dimension);

    Vec minB(std::numeric_limits<double>::max()), maxB(-std::numeric_limits<double>::max());
    for(int i=0; i<particleNumber; ++i)
    {
        if(particles[i].x[0]<minB[0]){minB[0]=particles[i].x[0];}
        if(particles[i].x[1]<minB[1]){minB[1]=particles[i].x[1];}
        if(particles[i].x[2]<minB[2]){minB[2]=particles[i].x[2];}
        if(particles[i].x[0]>maxB[0]){maxB[0]=particles[i].x[0];}
        if(particles[i].x[1]>maxB[1]){maxB[1]=particles[i].x[1];}
        if(particles[i].x[2]>maxB[2]){maxB[2]=particles[i].x[2];}
    }
    minB -= Vec(1.1*h);
    maxB += Vec(1.1*h);
    addBoundaryBox(minB, maxB, h);
    visuBox = Box(minB, maxB);

    //Neighbor Grid
    grid.minCorner = Vec( minB[0] - 20*h, minB[1] - 20*h, minB[2] - 20*h);
    grid.maxCorner = Vec( maxB[0] + 20*h, maxB[1] + 20*h, maxB[2] + 20*h);
    grid.init();
    for( int i = 0; i < boundaryNumber; ++i )
        grid.addBoundary(i, boundaries[i].x);
    for(int i = 0; i < particleNumber; ++i)
        grid.addFluid(i, particles[i]);
    mortonSort();

    //init
    computeBoundaryVolume();

    prepareGrid(); 

    for(int i=0; i<particleNumber; ++i)
    {
        computeRho(i);
        particles[i].rho = restDensity;
    }

    for(int i=0; i<particleNumber; ++i)
    {
        computeDii(i);
        //computeDii_Fluid(i);
        //computeDii_Boundary(i);
    }

    for(int i=0; i<particleNumber; ++i)
        computeAii(i);

    computeSceneParam();
    debugFluid();
}

void System::addBoundaryBox(const Vec& offset, const Vec& dimension)
{
    int epsilon = 0;
    int widthSize = floor(dimension[0]/h);
    int heightSize = floor(dimension[1]/h);
    int depthSize = floor(dimension[2]/h);

    //ZX plane - bottom
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec position(i*h, 0, j*h);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //ZX plane - top
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec position(i*h, dimension[1], j*h);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //XY plane - back
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= heightSize+epsilon; ++j)
        {
            Vec position(i*h, j*h, 0);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //XY plane - front
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= heightSize-epsilon; ++j)
        {
            Vec position(i*h, j*h, dimension[2]);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //YZ plane - left
    for(int i = -epsilon; i <= heightSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec position(0, i*h, j*h);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //YZ plane - right
    for(int i = -epsilon; i <= heightSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec position(dimension[0], i*h, j*h);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }
    translateBoundaries(offset);
    minBoundary = offset;
    maxBoundary = offset + dimension;
}

void System::createDamBreak()
{
    int wishedNumber = 1000;
    setParameters(wishedNumber, 1.50);

    p_kernel = MonaghanKernel( h );
    a_kernel = AkinciKernel( 2.0*h );
    b_kernel = BoundaryKernel( boundaryH, cs );

    Vec minOff(1.2*h,1.2*h,1.2*h);
    double size = pow(volume,1.0/3.0);
    createParticleVolume(minOff, 0.6*size, 1.0*size, 2.0*size, h, wishedNumber);

    Vec min(1000), max(0.0);
    for(int i=0; i<particleNumber; ++i)
    {
        if(particles[i].x[0]<min[0]){min[0]=particles[i].x[0];}
        if(particles[i].x[1]<min[1]){min[1]=particles[i].x[1];}
        if(particles[i].x[2]<min[2]){min[2]=particles[i].x[2];}
        if(particles[i].x[0]>max[0]){max[0]=particles[i].x[0];}
        if(particles[i].x[1]>max[1]){max[1]=particles[i].x[1];}
        if(particles[i].x[2]>max[2]){max[2]=particles[i].x[2];}
    }

    Vec minB = Vec(0,0,0);
    Vec maxB = Vec(2.0,2.0,2.35);
    addBoundaryBox(minB, maxB, h);
    visuBox = Box(minB, maxB);

    //Neighbor Grid
    grid.minCorner = Vec( minB[0] - 20*h, minB[1] - 20*h, minB[2] - 20*h);
    grid.maxCorner = Vec( maxB[0] + 20*h, maxB[1] + 20*h, maxB[2] + 20*h);
    grid.init();
    for( int i = 0; i < boundaryNumber; ++i )
        grid.addBoundary(i, boundaries[i].x);
    for(int i = 0; i < particleNumber; ++i)
        grid.addFluid(i, particles[i]);
    mortonSort();

    //init
    computeBoundaryVolume();

    prepareGrid(); 

    for(int i=0; i<particleNumber; ++i)
    {
        computeRho(i);
        particles[i].rho = restDensity;
    }

    for(int i=0; i<particleNumber; ++i)
    {
        computeDii(i);
        //computeDii_Fluid(i);
        //computeDii_Boundary(i);
    }

    for(int i=0; i<particleNumber; ++i)
        computeAii(i);

    computeSceneParam();
    debugFluid();
}


void System::translateBoundaries(const Vec& t)
{
    for(size_t i=0; i<boundaries.size(); ++i)
    {
        boundaries[i].x += t;
    }
}

void System::translateParticles(const Vec& t)
{
    for(size_t i=0; i<particles.size(); ++i)
    {
        particles[i].x += t;
    }
}

void System::addParticleMesh(const std::string& filename)
{
}

void System::addParticleBox(const Vec& offset, const Vec& dimension)
{
    int widthSize = floor(dimension[0]/h);
    int heightSize = floor(dimension[1]/h);
    int depthSize = floor(dimension[2]/h);

    for(int i=0; i<widthSize; ++i)
    {
        for(int j=0; j<heightSize; ++j)
        {
            for(int k=0; k<depthSize; ++k)
            {
                    Vec _x(i*h,j*h,k*h);
                    Vec _v(0,0,0);
                    particles.push_back( Particle(_x,_v) );
                    particleNumber++;
            }
        }
    }
    translateParticles(offset);
}

void System::addParticleBox(double width, double height, double depth, double spacing)
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
                    Vec _x(i*spacing,j*spacing,k*spacing);
                    Vec _v(0,0,0);
                    particles.push_back( Particle(_x,_v) );
                    particleNumber++;
            }
        }
    }
}


void System::createParticleVolume(Vec& pos, double width, double /*height*/, double depth, double spacing, int particleMax)
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
                    Vec _x(pos[0]+i*spacing,pos[1]+j*spacing,pos[2]+k*spacing);
                    Vec _v(0,0,0);
                    particles.push_back( Particle(_x,_v) );
                    particleNumber++;
                }
                count++;
            }
        }
        j++;
    }
}

void System::addBoundaryBox(const Vec& min, const Vec& max, const double spacing)
{
    double widthSize = floor( (max[0]-min[0])/(double)spacing );
    double heightSize = floor( (max[1]-min[1])/(double)spacing );
    double depthSize = floor( (max[2]-min[2])/(double)spacing );
    int epsilon = 0;

    //ZX plane - bottom
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec position(min[0] + i*spacing, min[1], min[2]+j*spacing);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //ZX plane - top
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec position(min[0] + i*spacing, max[1], min[2]+j*spacing);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //XY plane - back
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= heightSize+epsilon; ++j)
        {
            Vec position(min[0] + i*spacing, min[1]+j*spacing, min[2]);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //XY plane - front
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= heightSize-epsilon; ++j)
        {
            Vec position(min[0] + i*spacing, min[1]+j*spacing, max[2]);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //YZ plane - left
    for(int i = -epsilon; i <= heightSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec position(min[0], min[1]+i*spacing, min[2]+j*spacing);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //YZ plane - right
    for(int i = -epsilon; i <= heightSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec position(max[0], min[1]+i*spacing, min[2]+j*spacing);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }
}

bool pairCompare( const std::pair<int,int>& e1, const std::pair<int,int>& e2 )
{
    return (e1.second < e2.second);
}

void System::mortonSort()
{
    std::vector< std::pair<int, int> > particleZindex;

    //Fill particleZindex with particle index and Z-index
    for(int i = 0; i < particleNumber; ++i)
    {
        array<int,3> gridIndex = grid.getGridCoordinate( particles[i].x ); 
        int zindex = mortonNumber( gridIndex );
        std::pair<int,int> paire(i, zindex);
        particleZindex.push_back(paire);
    }

    //Sort according to z-index
    std::sort( particleZindex.begin(), particleZindex.end(), pairCompare );

    //Move particles according to z-index
    vector< Particle > oldParticles = particles;

    //double min = particleZindex[0].second;
    //double max = particleZindex[particleNumber-1].second;

    for(int i = 0; i < particleNumber; ++i)
    {
        std::pair<int,int>& paire = particleZindex[i];
        if(i != paire.first)
        {
            particles[i] = oldParticles[paire.first];
        }
        //double color = (paire.second-min)/(max-min);
        //particles[i].c[0] = color;
        //particles[i].c[1] = 0;
        //particles[i].c[2] = 0;
    }
}

void System::prepareGrid()
{
    if( countTime%100 == 0 )
        mortonSort();

    grid.clearFluid();
    for(int i = 0; i < particleNumber; ++i)
        grid.addFluid(i, particles[i]);

#pragma omp parallel for
    for(int i = 0; i < particleNumber; ++i)
        grid.getNearestNeighbor(i, particles, boundaries, 4.0*h*h);
}

void System::predictAdvection()
{
#pragma omp parallel for
    for(int i=0; i<particleNumber; ++i)
        computeRho(i);

#pragma omp parallel for
    for(int i=0; i<particleNumber; ++i)
        computeNormal(i);

#pragma omp parallel for
    for(int i=0; i<particleNumber; ++i)
    {
        computeAdvectionForces(i);
        predictVelocity(i);
        computeDii(i);
        //computeDii_Fluid(i);
        //computeDii_Boundary(i);
    }

#pragma omp parallel for
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
#pragma omp parallel for
        for(int i=0; i<particleNumber; ++i)
            computeSumDijPj(i);

#pragma omp parallel for
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

#pragma omp parallel for
    for(int i=0; i<particleNumber; ++i)
        computePressureForce(i);

#pragma omp parallel for
    for(int i=0; i<particleNumber; ++i)
    {
        Particle& pi=particles[i];
        pi.v = pi.v_adv + (dt*pi.f_p)/mass;
        pi.x += dt*pi.v;
    }
}

void System::simulate()
{
    prepareGrid();
    predictAdvection();
    pressureSolve();
    integration();
    computeStats();
}

//###################################################
//				System : Drawing function
//###################################################
void System::computeSceneParam()
{

    //Compute center
    scene_center = Vec(0,0,0);
    for(int i = 0; i < particleNumber; ++i)
    {
        scene_center += particles[i].x;
    }
    scene_center /= (double)particleNumber;

    //Compute radius
    scene_radius = 0.0f;
    Vec c;
    for(int i = 0; i < particleNumber; ++i)
    {
        c = particles[i].x-scene_center;
        double r = c.length();
        scene_radius = r>scene_radius ? r : scene_radius;
    }

    if(particleNumber < 10)
    {
        scene_center= Vec(0,0,0);
        scene_radius = 5;
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
    std::cout << "Cell Number : " << grid.cells.size() << std::endl;
    std::cout << "Size : " << grid.wNumber<< " x " << grid.hNumber<< " x " << grid.dNumber<< std::endl;
    std::cout << "Min : " << grid.minCorner[0] << ", " << grid.minCorner[1] << ", " << grid.minCorner[2] << std::endl;
    std::cout << "Max : " << grid.maxCorner[0] << ", " << grid.maxCorner[1] << ", " << grid.maxCorner[2] << std::endl;
}

/*
   void System::draw()
   {
//Draw grid
//grid.draw();
glDisable(GL_LIGHTING);
GLfloat black[] = {0.0f, 0.0f, 0.0f, 1.0f};
glMaterialfv(GL_FRONT, GL_DIFFUSE, black);
visuBox.draw();
glEnable(GL_LIGHTING);

//Draw fluid particle
//double radius = pow( 3.0*volume/(4*M_PI*particleNumber), 1.0/3.0 );
double radius = 0.5*h;
SolidSphere s(radius, 20,20);

//glDisable(GL_COLOR_MATERIAL);

//glColor3f(c[i][0], c[i][1], c[i][2]);
//GLfloat mat_diffuse[] = { 0.0f, 0.0f, 1.0f, 1.0f  };

GLfloat mat_ambient[] = { 1.0f, 1.0f, 1.0f, 1.0f };
GLfloat shininess[] = { 20.0f };
//GLfloat mat_emission[] =  { 0.3f, 0.2f, 0.2f, 0.0f };
GLfloat white[] = {1.0f, 1.0f, 1.0f, 1.0f};
glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
glMaterialfv(GL_FRONT, GL_SPECULAR, white);
glMaterialfv(GL_FRONT, GL_SHININESS, shininess);
//glMaterialfv(GL_FRONT, GL_EMISSION, mat_emission);
//glColor4f(0.0,0.0,1.0,1.0);

//glPointSize(5.0);
//glBegin(GL_POINTS);
for(int i = 0; i < particleNumber; ++i)
{
Particle& pi = particles[i];

GLfloat mat_diffuse[] = { pi.c[0], pi.c[1], pi.c[2], 1.0f  };
glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
//glVertex3f(pi.x[0], pi.x[1], pi.x[2]);
s.draw(pi.x[0], pi.x[1], pi.x[2]);
}
//glEnd();

//Draw boundary particles
//glColor4f(0.0,0.0,0.0,0.1f);
//glPointSize(0.1);
//GLfloat boundary_ambient[] = { 1.0f, 1.0f, 0.0f, 0.0f };
//GLfloat boundary_specular[] = { 1.0f, 1.0f, 1.0f, 0.0f };
//glMaterialfv(GL_FRONT, GL_AMBIENT, boundary_ambient);
//glMaterialfv(GL_FRONT, GL_SPECULAR, boundary_specular);
//for(int i = 0; i < boundaryNumber; ++i)
//{
//    //glVertex3f(boundaries[i].x[0], boundaries[i].x[1], boundaries[i].x[2]);
//    Boundary bi=boundaries[i];
//    s.draw(bi.x[0], bi.x[1], bi.x[2]);
//}
//glEnd();

glMaterialfv(GL_FRONT, GL_DIFFUSE, black);
}
*/

void System::write(const char * filename, vector<Vec3<double> > data)
{
    ofstream outputFile;
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
    vector< Vec3<double> > x = getPosition();
    vector< Vec3<double> > v = getVelocity();
    stringstream ss_padding; 
    ss_padding.fill('0');
    ss_padding.width(5);
    ss_padding << countExport++;
    std::string padding = ss_padding.str(); 

    stringstream posFilename, velFilename; 
    posFilename << baseName << "/position/position" << padding << ".txt";
    velFilename << baseName << "/velocity/velocity" << padding << ".txt";

    write(  posFilename.str().c_str(), x );
    write(  velFilename.str().c_str(), v );
}

