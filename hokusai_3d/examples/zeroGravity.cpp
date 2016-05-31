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

    FluidParams fluidParams = FluidParams(resolution, 1.0, 1000, 0.5, 1.0);
    BoundaryParams boundaryParams = BoundaryParams( fluidParams.smoothingRadius()/2.0, 0.0001, 1.0);

    Vec3r  fluidBox(1.5,1.5,1.5);
    Vec3r  fluidOffset(0,0,0);
    sph.addParticleBox(fluidOffset, fluidBox, fluidParams);

    Vec3r  boundBox(8.0,8.0,8.0);
    Vec3r  boundOffset = Vec3r(-4,-4,-4);
    sph.addBoundaryBox(boundOffset, boundBox, boundaryParams);

    Vec3r  gravity(0,0,0);
    sph.setGravity(gravity);

    sph.init();

    double time = 2.0;
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
