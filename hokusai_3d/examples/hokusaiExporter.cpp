#include <hokusai/system.hpp>
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
    int particleNumber = 1000; ///particle number
    HReal volume = 1.0; ///m3
    HReal restDensity = 1000.0; ///kg/m3
    HReal viscosity = 0.5;
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

    Vec3r  fluidBox(1.5,1.5,1.5);
    Vec3r  fluidOffset(0,0,0);
    Vec3r velocity(0,0,0);
    sph.addParticleBox(fluidOffset, fluidBox, velocity);

    Vec3r  gravity(0,0,0);
    sph.setGravity(gravity);

    sph.init();

    double time = 1.0;
    boost::timer::auto_cpu_timer t;
    boost::progress_display show_progress( std::floor(time/solverParams.timeStep()) );
    while(sph.getTime()<=time)
    {
        //Simulate
        sph.computeSimulationStep();
        //Update progress bar
        ++show_progress;
    }
    std::string fileName = "./../../data/hokusaiExporter.bin";
    HokusaiExporter hokusaiExporter(fileName, sph);

    return 0;
}
