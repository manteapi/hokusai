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
    int resolution = 10000; ///particle number per m3
    System sph(resolution);

    Vec fluidBox(1.0,2.0,1.0);
    Vec fluidOffset(0,0,0);
    sph.addParticleBox(fluidOffset, fluidBox);

    Vec securityOffset(1.05*sph.getSmoothingRadius());
    Vec boundBox(2.5,2.5,1.0);
    boundBox += securityOffset;
    Vec boundOffset = fluidOffset;
    boundOffset -= securityOffset;
    sph.addBoundaryBox(boundOffset, boundBox);

    sph.init();

    sph.getViscosity() = 0.1;
    sph.getFluidCohesion() = 0.05;
    sph.getBoundaryAdhesion() = 0.0001;

    std::cout << "Viscosity : " << sph.getViscosity() << std::endl;
    std::cout << "Cohesion : " << sph.getFluidCohesion() << std::endl;
    std::cout << "Adhesion : " << sph.getBoundaryAdhesion() << std::endl;

    double time = 5.0;
    int count=0;
    boost::timer::auto_cpu_timer t;
    boost::progress_display show_progress( std::floor(time/sph.getTimeStep()) );
    while(sph.getTime()<=time)
    {
        //Simulate
        sph.simulate();

        //Output
        if( std::floor((sph.getTime()-sph.getTimeStep())/0.016) != std::floor(sph.getTime()/0.016) )
        {
            write_frame(sph.particles, count, 10.0);
            sph.exportState("./");
            ++count;
        }

        //Update progress bar
        ++show_progress;
    }

    return 0;
}
