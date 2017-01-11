#include <hokusai/system.hpp>
#include <hokusai/io.hpp>
#include <hokusai/particleSource.hpp>

#define timer   timer_class
#include <boost/progress.hpp>
#undef timer
#include <boost/timer/timer.hpp>

#include <iostream>
#include <algorithm>

using namespace std;
using namespace hokusai;

int main()
{
    int particleNumber = 1e5; ///particle number
    HReal volume = 1.0; ///m3
    HReal restDensity = 1000.0; ///kg/m3
    HReal viscosity = 1e-5;
    HReal cohesion = 1e-3;
    FluidParams fluidParams(particleNumber, volume, restDensity, viscosity, cohesion);

    HReal adhesion=0.000;
    HReal friction=0.00;
    HReal samplingResolution = 0.2;
    BoundaryParams boundaryParams(samplingResolution*fluidParams.smoothingRadius(), adhesion, friction);

    HReal timeStep = 5e-4;
    int maxPressureSolveIterationNb = 2;
    HReal maxDensityError = 1.0;
    SolverParams solverParams(timeStep, maxPressureSolveIterationNb, maxDensityError);

    System sph(fluidParams, boundaryParams, solverParams);

    std::string meshFilename = "./../../mesh/suzanne.obj";
    sph.addBoundaryMesh(meshFilename.c_str());

    /*
       Vec3r offsetSphere(0,0,0);
       double radius = 0.3;
       Vec3r velocity(0,0,0);
       sph.addParticleSphere(offsetSphere, radius, velocity);
       */

    HReal startTime, endTime, delay, spacing;
    Vec3r  position, orientation, scale, velocity;

    startTime=0;
    endTime = 5.0;
    delay=0.010;
    spacing = 1.05*fluidParams.smoothingRadius();
    velocity = Vec3r(0.0,0.0,2.4);
    scale = Vec3r(0.2,0.2,0.2);

    position = Vec3r(0.0,0.8,0.0);
    orientation = Vec3r(M_PI/2.0,0,0);
    ParticleSource source(startTime, endTime, delay, spacing, position, orientation, scale, velocity);
    sph.addParticleSource(source);


    sph.init();

    double time = 7;
    int count=0;
    int frameNumber = std::floor(time/0.016)-1;
    boost::timer::auto_cpu_timer t;
    boost::progress_display show_progress( std::floor(time/solverParams.timeStep()) );
    std::string prefix="pouring_Ferstl2016";
    std::string path="./";
    while(sph.getTime()<=time)
    {
        //Simulate
        sph.computeSimulationStep();

        //Output
        if( std::floor((sph.getTime()-solverParams.timeStep())/0.016) != std::floor(sph.getTime()/0.016) )
        {
            write_frame(sph, count, 10.0);
            ++count;
        }

        //Update progress bar
        ++show_progress;
    }

    return 0;
}
