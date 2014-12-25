#include "../include/hokusai/system.hpp"
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

    gravity = Vec(0,-9.81,0);

    a_kernel = AkinciKernel();
    p_kernel = MonaghanKernel();
    b_kernel = BoundaryKernel();

    gridInfo = GridUtility();
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

    gravity = Vec(0,-9.81,0);

    a_kernel = AkinciKernel();
    p_kernel = MonaghanKernel();
    b_kernel = BoundaryKernel();

    gridInfo = GridUtility();
    particles = vector<Particle>();
    boundaries = vector<Boundary>();

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

void System::getNearestNeighbor(const int i, const float radius)
{
    Particle& p = particles[i];
    p.fluidNeighbor.clear();
    p.boundaryNeighbor.clear();

    std::vector<int> neighborCell;
    gridInfo.get27Neighbors(neighborCell, p.x, radius);

    for(size_t i=0; i<neighborCell.size(); ++i)
    {
        std::vector<int>& bNeighborCell = boundaryGrid[neighborCell[i]];
        for(size_t j=0; j<bNeighborCell.size(); ++j)
        {
            int bParticleId = boundaryGrid[neighborCell[i]][j];
            Boundary& bParticle = boundaries[bParticleId];
            Vec d = bParticle.x-p.x;
            if( d.lengthSquared()<radius*radius )
                p.boundaryNeighbor.push_back(bParticleId);
        }

        std::vector<int>& fNeighborCell = fluidGrid[neighborCell[i]];
        for(size_t j=0; j< fNeighborCell.size(); ++j)
        {
            int fParticleId = fluidGrid[neighborCell[i]][j];
            Particle& fParticle = particles[fParticleId];
            Vec d = fParticle.x-p.x;
            if( d.lengthSquared()<radius*radius)
                p.fluidNeighbor.push_back(fParticleId);
        }
    }
}

void System::getNearestNeighbor(vector< int >& neighbor, const vector< vector<int> >& grid, const Vec &x)
{
    std::vector<int> neighborCell;
    gridInfo.get27Neighbors(neighborCell, x, gridInfo.spacing());
    neighbor.clear();
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
    for(int i=0; i<boundaryNumber; ++i)
    {
        double densityNumber=0.0;
        vector<int> neighbors;
        getNearestNeighbor(neighbors, boundaryGrid, boundaries[i].x);
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
    boundaryH = h/2.0; //boundaryH must be <= h (neighbor search purpose)

    dt = 0.004;

    p_kernel = MonaghanKernel( h );
    a_kernel = AkinciKernel( 2.0*h );
    b_kernel = BoundaryKernel( boundaryH, cs );
}

void System::init()
{  
    mortonSort();

    boundaryGrid.resize(gridInfo.size());
    for(size_t i=0; i<boundaries.size(); ++i)
    {
        int id = gridInfo.cellId(boundaries[i].x);
        if(gridInfo.isInside(id))
            boundaryGrid[id].push_back(i);
    }

    fluidGrid.resize(gridInfo.size());
    for(size_t i=0; i<particles.size(); ++i)
    {
        int id = gridInfo.cellId(particles[i].x);
        if(gridInfo.isInside(id))
            fluidGrid[id].push_back(i);
    }

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

    debugFluid();
}

void System::addBoundaryBox(const Vec& offset, const Vec& scale)
{
    int epsilon = 0;
    int widthSize = floor(scale[0]/h);
    int heightSize = floor(scale[1]/h);
    int depthSize = floor(scale[2]/h);

    //ZX plane - bottom
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec position(i*h, offset[1], j*h);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //ZX plane - top
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec position(i*h, offset[1]+scale[1], j*h);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //XY plane - back
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= heightSize+epsilon; ++j)
        {
            Vec position(i*h, j*h, offset[2]);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //XY plane - front
    for(int i = -epsilon; i <= widthSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= heightSize-epsilon; ++j)
        {
            Vec position(i*h, j*h, offset[2]+scale[2]);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //YZ plane - left
    for(int i = -epsilon; i <= heightSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec position(offset[0], i*h, j*h);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    //YZ plane - right
    for(int i = -epsilon; i <= heightSize+epsilon; ++i)
    {
        for(int j = -epsilon; j <= depthSize+epsilon; ++j)
        {
            Vec position(offset[0]+scale[0], i*h, j*h);
            boundaries.push_back(Boundary(position,Vec(0.0),0.0));
            boundaryNumber++;
        }
    }

    gridInfo.update(offset-Vec(2.0*h), scale+Vec(4.0*h), 2.0*h);
}

void System::addBoundarySphere(const Vec& offset, const double& radius)
{
    std::vector<Vec> samples = getSphereSampling(offset, radius, h, h);
    for(size_t i=0; i<samples.size(); ++i)
    {
        boundaries.push_back(Boundary(samples[i],Vec(0.0),0.0));
        boundaryNumber++;
    }
}

void System::addBoundaryHemiSphere(const Vec& offset, const double& radius)
{
    std::vector<Vec> samples = getHemiSphereSampling(offset, radius, h, h);
    for(size_t i=0; i<samples.size(); ++i)
    {
        boundaries.push_back(Boundary(samples[i],Vec(0.0),0.0));
        boundaryNumber++;
    }
}

void System::addBoundaryDisk(const Vec& offset, const double& radius)
{
    std::vector<Vec> samples = getDiskSampling(offset, radius, h);
    for(size_t i=0; i<samples.size(); ++i)
    {
        boundaries.push_back(Boundary(samples[i],Vec(0.0),0.0));
        boundaryNumber++;
    }
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

void System::addParticleSphere(const Vec& centre, const double radius)
{
    Vec scale(2.0*radius, 2.0*radius, 2.0*radius);
    Vec offset = centre - Vec(radius, radius, radius);
    GridUtility grid(offset, scale, h);

    for(int i=0; i<grid.size(); ++i)
    {
        Vec _x = grid.gridToWorld(i);
        _x+=grid.spacing()/2.0;
        Vec _v(0,0,0);
        double l2 = (centre-_x).lengthSquared();
        if(l2<=(radius*radius))
        {
            particles.push_back( Particle(_x,_v) );
            particleNumber++;
        }
    }
}

void System::addParticleBox(const Vec& offset, const Vec& scale)
{
    int widthSize = floor(scale[0]/h);
    int heightSize = floor(scale[1]/h);
    int depthSize = floor(scale[2]/h);

    for(int i=0; i<widthSize; ++i)
    {
        for(int j=0; j<heightSize; ++j)
        {
            for(int k=0; k<depthSize; ++k)
            {
                Vec _x = offset + Vec(i*h,j*h,k*h);
                Vec _v(0,0,0);
                particles.push_back( Particle(_x,_v) );
                particleNumber++;
            }
        }
    }
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

void System::addFluidParticle(const Vec& x, const Vec& v)
{
    particles.push_back( Particle(x,v) );
    particleNumber++;
}

void System::addBoundaryParticle(const Vec& x, const Vec& v)
{
    boundaries.push_back( Boundary(x,v) );
    boundaryNumber++;
}

void System::addBoundaryMesh(const char* filename)
{
    TriMesh mesh(filename);
    std::vector<Vec3f> samples;
    AkinciMeshSampling(mesh, h/2.0, samples);
    Vec3f minBB(std::numeric_limits<double>::max()), maxBB(-std::numeric_limits<double>::max());
    for(size_t i=0; i<samples.size(); ++i)
    {
        for(int j=0; j<3; ++j)
        {
            minBB[j] = std::min(samples[i][j], minBB[j]);
            maxBB[j] = std::max(samples[i][j], maxBB[j]);
        }
        boundaries.push_back(Boundary(samples[i],Vec(0.0),0.0));
        boundaryNumber++;
    }
    Vec3f offset = minBB;
    Vec3f scale = maxBB-minBB;
    gridInfo.update(offset-Vec(2.0*h), scale+Vec(4.0*h), 2.0*h);
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
        Vec3i _gridIndex = gridInfo.worldToGrid(particles[i].x);
        array<int,3> gridIndex;
        for(int j=0; j<3; ++j)
            gridIndex[j] = _gridIndex[j];
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

    for(size_t i=0; i<fluidGrid.size(); ++i)
        fluidGrid[i].clear();

    for(size_t i=0; i<particles.size(); ++i)
    {
        int id = gridInfo.cellId(particles[i].x);
        if(gridInfo.isInside(id))
            fluidGrid[id].push_back(i);
    }

#ifdef HOKUSAI_USING_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < particleNumber; ++i)
        getNearestNeighbor(i, 2.0*h);
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

void System::simulate()
{
    prepareGrid();
    predictAdvection();
    pressureSolve();
    integration();
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
    gridInfo.info();
}

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

}
