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

    int resolution = 100; ///particle number per m3
    System sph(resolution);

    FluidParams fluidParams = FluidParams(resolution, 1.0, 1000, 0.1, 0.05);
    BoundaryParams boundaryParams = BoundaryParams( fluidParams.smoothingRadius()/2.0, 0.0001, 1.0);

    Vec3r  offsetDomain(0,0,0);
    Vec3r  scaleDomain(6,6,6);
    sph.addBoundaryBox(offsetDomain, scaleDomain, boundaryParams);

    Vec3r  offsetSphere(3,3,3);
    double radius = 2.5;
    sph.addParticleSphere(offsetSphere, radius, fluidParams);

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
