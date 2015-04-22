#include <hokusai/system.hpp>
#include <hokusai/utils.hpp>
#include <hokusai/particleSource.hpp>

#define timer   timer_class
#include <boost/progress.hpp>
#undef timer
#include <boost/timer/timer.hpp>

#include <iostream>
#include <algorithm>

using namespace std;
using namespace hokusai;

typedef double SReal;

int main()
{

    int resolution = 1000; ///particle number per m3
    System sph(resolution);

//    Vec fluidBox(1.0,2.0,1.0);
//    Vec fluidOffset(0,0,0);
//    sph.addParticleBox(fluidOffset, fluidBox);

    Vec securityOffset(1.05*sph.getSmoothingRadius());
    Vec boundBox(2.5,2.5,1.0);
    boundBox += securityOffset;
    Vec boundOffset(0,0,0);
    boundOffset -= securityOffset;
    sph.addBoundaryBox(boundOffset, boundBox);

    SReal startTime, endTime, delay;
    Vec position, orientation, scale, velocity;

    startTime=0;
    endTime = 10;
    delay=0.04;
    position = Vec(0.5,0.5,0.5);
    orientation = Vec(0,0,0);
    scale = Vec(1,1,1);
    velocity = Vec(0.1,0,0);

    ParticleSource source(startTime, endTime, delay, position, orientation, scale, velocity);
    sph.addParticleSource(source);

    sph.init();

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
            write_frame(sph.particles, count);
            sph.exportState("./");
            ++count;
        }

        //Update progress bar
        ++show_progress;
    }
    std::cout << "Particle Number : " << sph.getParticleNumber() << std::endl;

    return 0;
}
