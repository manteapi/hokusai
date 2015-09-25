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
    int resolution = 20000; ///particle number per m3
    System sph(resolution);

    Vec securityOffset(1.05*sph.getSmoothingRadius());
    Vec boundBox(2.0,2.5,1.0);
    boundBox += securityOffset;
    Vec boundOffset(0,0,0);
    boundOffset -= securityOffset;
    sph.addBoundaryBox(boundOffset, boundBox);

    SReal startTime, endTime, delay, spacing;
    Vec position, orientation, scale, velocity;

    startTime=0;
    endTime = 3.0;
    delay=0.02;
    spacing = 1.05*sph.getSmoothingRadius();
    velocity = Vec(0.0,0.0,2.0);
    scale = Vec(0.15,0.15,0.15);

    position = Vec(0.25,1.0,0.5);
    orientation = Vec(0,M_PI/2.0,0);
    ParticleSource source1(startTime, endTime, delay, spacing, position, orientation, scale, velocity);
    sph.addParticleSource(source1);

    position = Vec(1.75,1.0,0.5);
    orientation = Vec(0,-M_PI/2.0,0);
    ParticleSource source2(startTime, endTime, delay, spacing, position, orientation, scale, velocity);
    sph.addParticleSource(source2);

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
