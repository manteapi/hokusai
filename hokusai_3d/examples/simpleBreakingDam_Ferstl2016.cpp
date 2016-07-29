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
    int particleNumber = 1.6e6; ///particle number
    HReal volume = 1.0; ///m3
    HReal restDensity = 1000.0; ///kg/m3
    HReal viscosity = 1e-5;
    HReal cohesion = 1e-3;
    FluidParams fluidParams(particleNumber, volume, restDensity, viscosity, cohesion);

    HReal adhesion=0.000;
    HReal friction=0.00;
    BoundaryParams boundaryParams(0.4*fluidParams.smoothingRadius(), adhesion, friction);

    HReal timeStep = 1e-4;
    int maxPressureSolveIterationNb = 2;
    HReal maxDensityError = 1.0;
    SolverParams solverParams(timeStep, maxPressureSolveIterationNb, maxDensityError);

    System sph(fluidParams, boundaryParams, solverParams);

    Vec3r securityOffset(1.05*fluidParams.smoothingRadius());
    Vec3r boundaryOffset(0.0,0.0,0.0);
    Vec3r boundBox(1.0,1.0,1.0);
    boundBox += 2.0*securityOffset;
    boundaryOffset -= securityOffset;
    sph.addBoundaryBox(boundaryOffset, boundBox);

    //Import
    /*
    std::string fileName = "./../../data/simpleBreakingDam_Ferstl2016.bin";
    HokusaiImporter hokusaiImporter(fileName, sph);
    */

    //Export
    Vec3r boundaryOffset2(0.6,0.2,0.0);
    Vec3r boundBox2(0.4,0.8,1.0);
    sph.addBoundaryBox(boundaryOffset2, boundBox2);

    Vec3r boundaryOffset3(0.0,0.2,0.6);
    Vec3r boundBox3(1.0,0.8,0.4);
    sph.addBoundaryBox(boundaryOffset3, boundBox3);

    Vec3r fluidBox1(0.6,0.6,0.6);
    fluidBox1 = fluidBox1;
    Vec3r fluidOffset1(0,0,0);
    Vec3r velocity1(0,0,0);
    sph.addParticleBox(fluidOffset1, fluidBox1, velocity1);

    Vec3r fluidBox2(0.4,0.2,1.0);
    fluidBox2 = fluidBox2;
    Vec3r fluidOffset2(0.6+0.0*fluidParams.smoothingRadius(),0,0);
    Vec3r velocity2(0,0,0);
    sph.addParticleBox(fluidOffset2, fluidBox2, velocity2);

    Vec3r fluidBox3(0.6,0.2,0.4);
    Vec3r fluidOffset3(0.0,0.0,0.6+0.0*fluidParams.smoothingRadius());
    Vec3r velocity3(0,0,0);
    sph.addParticleBox(fluidOffset3, fluidBox3, velocity3);

    sph.init();

    double time = 10;
    int count=0;
    HReal framerate = 1/25.0;
    int frameNumber = std::floor(time/framerate)-1;
    boost::timer::auto_cpu_timer t;
    boost::progress_display show_progress( std::floor(time/solverParams.timeStep()) );
    std::string prefix="/media/manteapi/49534769-eb68-4ca6-a09a-7e173f850b03/Hokusai/simpleBreakingDam_Ferstl2016_2x/simpleBreakingDam_Ferstl2016";
    std::string path="./";
    BlenderExporter blenderExporter(prefix, path, sph.particleNumber(), frameNumber);
    while(sph.getTime()<=time)
    {
        //Simulate
        sph.computeSimulationStep();

        //Output
        if( std::floor((sph.getTime()-solverParams.timeStep())/framerate) != std::floor(sph.getTime()/framerate) )
        {
            blenderExporter.apply(sph);
            write_frame(sph, count, 10.0);
            ++count;
        }

        //Update progress bar
        ++show_progress;
    }

    std::string outputFileName= "./../../data/simpleBreakingDam_Ferstl2016_2x.bin";
    //std::string outputFileName = "./../../data/save.bin";
    HokusaiExporter hokusaiExporter(outputFileName, sph);

    return 0;
}
