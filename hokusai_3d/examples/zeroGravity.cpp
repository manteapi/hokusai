#include <hokusai/solver/iisphSolver.inl>
#include <hokusai/system.inl>
#include <hokusai/utils.inl>

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
    System<IISPHSolver> sph(resolution);

    Vec3r  fluidBox(1.5,1.5,1.5);
    Vec3r  fluidOffset(0,0,0);
    sph.addParticleBox(fluidOffset, fluidBox);

    Vec3r  boundBox(8.0,8.0,8.0);
    Vec3r  boundOffset = Vec3r(-4,-4,-4);
    sph.addBoundaryBox(boundOffset, boundBox);

    Vec3r  gravity(0,0,0);
    sph.setGravity(gravity);

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
            write_frame< System<IISPHSolver>::Particle >(sph.getParticles(), count);
            ++count;
        }

        //Update progress bar
        ++show_progress;
    }

    return 0;
}
