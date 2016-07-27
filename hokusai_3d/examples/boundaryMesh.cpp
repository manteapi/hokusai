#include <hokusai/system.hpp>
#include <hokusai/fluidParams.hpp>
#include <hokusai/io.hpp>

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
    int particleNumber = 2000; ///particle number
    HReal volume = 1.0; ///m3
    HReal restDensity = 1000.0; ///kg/m3
    HReal viscosity = 0.1;
    HReal cohesion = 0.05;
    FluidParams fluidParams(particleNumber, volume, restDensity, viscosity, cohesion);

    HReal adhesion=0.001;
    HReal friction=1.0;
    BoundaryParams boundaryParams(fluidParams.smoothingRadius()/2.0, adhesion, friction);

    HReal timeStep = 4e-3;
    int maxPressureSolveIterationNb = 2;
    HReal maxDensityError = 1.0;
    SolverParams solverParams(timeStep, maxPressureSolveIterationNb, maxDensityError);

    System sph(fluidParams, boundaryParams, solverParams);

    std::string filename = "./../../mesh/sphere.obj";
    sph.addBoundaryMesh(filename.c_str());

    sph.m_gridInfo.info();

    Vec3r offsetSphere(1,1,0);
    double radius = 0.5;
    Vec3r velocity(0,0,0);
    sph.addParticleSphere(offsetSphere, radius, velocity);

    sph.init();

    double time = 2.0;
    int count=0;
    boost::timer::auto_cpu_timer t;
    boost::progress_display show_progress( std::floor(time/solverParams.timeStep()) );
    while(sph.getTime()<=time)
    {
        //Simulate
        sph.computeSimulationStep();

        //Output
        if( std::floor((sph.getTime()-solverParams.timeStep())/0.016) != std::floor(sph.getTime()/0.016) )
        {
            sph.exportState("./output/");
            write_frame(sph, count);
            ++count;
        }

        //Update progress bar
        ++show_progress;
    }

    return 0;
}
