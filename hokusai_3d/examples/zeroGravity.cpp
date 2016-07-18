#include <hokusai/system.hpp>
#include <hokusai/utils.hpp>
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
    int particleNumber = 4e4; ///particle number
    HReal volume = 1.0; ///m3
    HReal restDensity = 1000.0; ///kg/m3
    HReal viscosity = 1e-5;
    HReal cohesion = 0.05;
    FluidParams fluidParams(particleNumber, volume, restDensity, viscosity, cohesion);
    HReal adhesion=0.001;
    HReal friction=1.0;
    BoundaryParams boundaryParams(fluidParams.smoothingRadius()/2.0, adhesion, friction);

    HReal timeStep = 1e-2;
    int maxPressureSolveIterationNb = 2;
    HReal maxDensityError = 1.0;
    SolverParams solverParams(timeStep, maxPressureSolveIterationNb, maxDensityError);

    System sph(fluidParams, boundaryParams, solverParams);

    HReal radius1 = 0.5;
    Vec3r  fluidOffset1(0,0,0);
    Vec3r velocity1(0.1,0,0);
    sph.addParticleSphere(fluidOffset1, radius1, velocity1);

    HReal radius2 = 0.5;
    Vec3r  fluidOffset2(1.5,0,0);
    Vec3r velocity2(-0.1,0,0);
    sph.addParticleSphere(fluidOffset2, radius2, velocity2);

    Vec3r  gravity(0,0,0);
    sph.setGravity(gravity);

    sph.init();

    double time = 240.0;
    int count=0;
    int frameNumber = std::floor(time/0.016)-1;

    std::string prefix="test";
    std::string path="/media/manteapi/DATA/FluidSimulation/";
    BlenderExporter blenderExporter(prefix, path, sph.particleNumber(), frameNumber);

    boost::timer::auto_cpu_timer t;
    boost::progress_display show_progress( std::floor(time/solverParams.timeStep()) );
    while(sph.getTime()<=time)
    {
        //Simulate
        sph.computeSimulationStep();

        //Output
        if( std::floor((sph.getTime()-solverParams.timeStep())/0.016) != std::floor(sph.getTime()/0.016) )
        {
            blenderExporter.apply(sph);
            write_frame(sph.m_particles, count);
            ++count;
        }

        //Update progress bar
        ++show_progress;
    }

    return 0;
}

    HReal timeStep = 1e-2;
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

    double time = 60.0;
    int count=0;
    int frameNumber = std::floor(time/0.016)-1;

    std::string prefix="test";
    std::string path="/media/manteapi/49534769-eb68-4ca6-a09a-7e173f850b03/Hokusai/ZeroGravity/"
    BlenderExporter blenderExporter(prefix, path, sph.particleNumber(), frameNumber);

    boost::timer::auto_cpu_timer t;
    boost::progress_display show_progress( std::floor(time/solverParams.timeStep()) );
    while(sph.getTime()<=time)
    {
        //Simulate
        sph.computeSimulationStep();

        //Output
        if( std::floor((sph.getTime()-solverParams.timeStep())/0.016) != std::floor(sph.getTime()/0.016) )
        {
            blenderExporter.apply(sph);
            write_frame(sph.m_particles, count);
            ++count;
        }

        //Update progress bar
        ++show_progress;
    }

    return 0;
}
