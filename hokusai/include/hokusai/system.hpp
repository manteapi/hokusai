#ifndef HOKUSAI_SYSTEM_HPP
#define HOKUSAI_SYSTEM_HPP

#include <iostream>
#include <fstream>
#include <ctime>

#include "Vec.hpp"
#include "utils.hpp"
#include "particle.hpp"
#include "gridUtility.hpp"
#include "triMesh.hpp"
#include "sampler.hpp"


namespace hokusai
{

typedef Vec3<double> Vec;
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

    double volume;
    double restDensity;
    double mean_density;
    double density_fluctuation;
    double real_volume;
    double mass;
    double h; // Smoothing radius
    double fcohesion;
    double badhesion;
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

public :
    void getNearestNeighbor(vector< int >& neighbors, const vector<vector<int> > &grid, const Vec& x);
    void getNearestNeighbor(const int i, const float radius);

    //Simulation Loop
    void prepareGrid();
    void predictAdvection();
    void predictRho(int i);
    void initializePressure(int i);
    void computeNormal(int i);
    bool isSurfaceParticle(int i, double treshold);
    vector<Particle> getSurfaceParticle();
    void computeRho(int i);
    void computeAdvectionForces(int i);
    void computeViscosityForces(int i, int j);
    void computeSurfaceTensionForces(int i, int j);
    void computeBoundaryFrictionForces(int i, int j);
    void computeBoundaryAdhesionForces(int i, int j);
    void predictVelocity(int i);
    void computeDii(int i);
    void computeDii_Fluid(int i);
    void computeDii_Boundary(int i);
    void computeAii(int i);
    void pressureSolve();
    void computeSumDijPj(int i);
    Vec computeDij(int i, int j);
    void computePressure(int i);
    void computePressureForce(int i);
    void computeError();
    void integration();
    void simulate();

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

    void addBoundaryMesh(const char* filename);
    void addBoundaryBox(const Vec& min, const Vec& max, const double spacing);
    void addBoundaryBox(const Vec& min, const Vec& scale);

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
    vector< Vec3<double> > getPosition(){ vector<Vec3<double> > pos; for(int i=0; i<particleNumber; ++i){pos.push_back(particles[i].x);} return pos;}
    vector< Vec3<double> > getVelocity(){ vector<Vec3<double> > vel; for(int i=0; i<particleNumber; ++i){vel.push_back(particles[i].v);} return vel;}
    vector< Vec3<double> > getNormal(){ vector<Vec3<double> > normal; for(int i=0; i<particleNumber; ++i){normal.push_back(particles[i].n);} return normal;}
    void write(const char* filename, vector< Vec3<double> > data);
    void exportState(const char* baseName);
    double getTimeStep(){return dt;}
    void setTimeStep(const double _dt){dt = _dt;}
    double getTime(){return time;}
    double & getSmoothingRadiusValue(){ return h; }
    const double & getSmoothingRadius(){ return h; }
    double & getTimeStepValue(){ return dt; }
    double & getMassValue(){ return mass; }
    double & getMeanDensityValue(){ return mean_density;}
    double & getDensityFluctuationValue(){ return density_fluctuation;}
    double & getRealVolumeValue(){ return real_volume; }
    int & getParticleNumber(){ return particleNumber; }

    void applyShepardFilter();
};

bool pairCompare( const std::pair<int,int>& e1, const std::pair<int,int>& e2 );
}

#endif // SYSTEM_H
