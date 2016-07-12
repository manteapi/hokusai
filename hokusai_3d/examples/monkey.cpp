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

    int resolution = 1000; ///particle number per m3
    System sph(resolution);

    Vec3r  fluidBox(1.0,2.0,1.0);
    Vec3r  fluidOffset(0,0,0);
    sph.addParticleBox(fluidOffset, fluidBox);

    Vec3r  securityOffset(1.05*sph.getSmoothingRadius());
    Vec3r  boundBox(2.5,2.5,1.0);
    boundBox += securityOffset;
    Vec3r  boundOffset = fluidOffset;
    boundOffset -= securityOffset;
    sph.addBoundaryBox(boundOffset, boundBox);

    sph.init();

    double time = 1.0;
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
            write_frame(sph.m_particles, count);
            ++count;
        }

        //Update progress bar
        ++show_progress;
    }

    return 0;
}
