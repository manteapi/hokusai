#include <hokusai/solver/iisphSolver.inl>
#include <hokusai/system.inl>
#include <hokusai/utils.inl>
#include <hokusai/grid.hpp>

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

    int resolution = 2000; ///particle number per m3
    System sph(resolution);

    Vec3r  fluidBox(2.0,4.0,1.0);
    Vec3r  fluidOffset(0,0,0);
    sph.addParticleBox(fluidOffset, fluidBox);

    Vec3r  securityOffset(1.05*sph.getSmoothingRadius());
    Vec3r  boundBox(6.0,8.0,1.0);
    boundBox += securityOffset;
    Vec3r  boundOffset = fluidOffset;
    boundOffset -= securityOffset;
    sph.addBoundaryBox(boundOffset, boundBox);

    sph.getTimeStep() = 0.002;

    sph.init();

    double time = 6;
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
            //write_frame< ParticleIISPH >(sph.getParticles(), count);
            sph.exportState("./output");
            ++count;
        }

        //Update progress bar
        ++show_progress;
    }

    return 0;
}
