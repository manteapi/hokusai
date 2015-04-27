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
    int resolution = 290000; ///particle number per m3
    System sph(resolution);

    Vec fluidBox1(0.5,1.0,0.5);
    Vec fluidOffset1(0,0,0);
    sph.addParticleBox(fluidOffset1, fluidBox1);

    Vec fluidBox2(0.5,1.0,0.5);
    Vec fluidOffset2(1.0,0.0,0.5);
    sph.addParticleBox(fluidOffset2, fluidBox2);

    Vec securityOffset(1.05*sph.getSmoothingRadius());
    Vec boundBox(1.5,1.0,1.0);
    boundBox += securityOffset;
    Vec boundOffset = fluidOffset1;
    boundOffset -= securityOffset;
    sph.addBoundaryBox(boundOffset, boundBox);

    sph.init();

    sph.getViscosity() = 0.01;
    sph.getFluidCohesion() = 1.0;
    sph.getBoundaryAdhesion() = 0.01;
    sph.getBoundaryFriction() = 0.01;
    sph.getTimeStep() = 1e-3;

    std::cout << "Viscosity : " << sph.getViscosity() << std::endl;
    std::cout << "Cohesion : " << sph.getFluidCohesion() << std::endl;
    std::cout << "Adhesion : " << sph.getBoundaryAdhesion() << std::endl;

    double time = 10.0;
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
