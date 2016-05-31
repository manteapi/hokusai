#include <hokusai/system.hpp>
#include <hokusai/utils.hpp>

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
    int resolution = 2.5e6; ///particle number per m3
    System sph(resolution);

    FluidParams fluidParams = FluidParams(resolution, 1.0, 1000, 0.1, 0.05);
    BoundaryParams boundaryParams = BoundaryParams( fluidParams.smoothingRadius()/2.0, 0.0001, 1.0);

    Vec3r  fluidBox1(0.5,1.0,0.5);
    Vec3r  fluidOffset1(0,0,0);
    sph.addParticleBox(fluidOffset1, fluidBox1, fluidParams);

    Vec3r  fluidBox2(0.5,1.0,0.5);
    Vec3r  fluidOffset2(1.0,0.0,0.5);
    sph.addParticleBox(fluidOffset2, fluidBox2, fluidParams);

    Vec3r  securityOffset(1.05*fluidParams.smoothingRadius());
    Vec3r  boundBox(1.5,1.0,1.0);
    boundBox += securityOffset;
    Vec3r  boundOffset = fluidOffset1;
    boundOffset -= securityOffset;
    sph.addBoundaryBox(boundOffset, boundBox,boundaryParams);

    sph.init();
    sph.getTimeStep() = 1e-4;

    double time = 10.0;
    int count=0;
    boost::timer::auto_cpu_timer t;
    boost::progress_display show_progress( std::floor(time/sph.getTimeStep()) );
    while(sph.getTime()<=time)
    {
        //Simulate
        sph.computeSimulationStep();

        //Output
        if( std::floor((sph.getTime()-sph.getTimeStep())/0.016) != std::floor(sph.getTime()/0.016) )
        {
            write_frame(sph.m_particles, count, 10.0);
            sph.exportState("./");
            ++count;
        }

        //Update progress bar
        ++show_progress;
    }

    return 0;
}
