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

    int resolution = 2000; ///particle number per m3
    System sph(resolution);

    std::string filename = "./../mesh/sphere.obj";
    sph.addBoundaryMesh(filename.c_str());

    sph.gridInfo.info();

    Vec offsetSphere(1,1,0);
    double radius = 0.5;
    sph.addParticleSphere(offsetSphere, radius);

    sph.init();

    double time = 2.0;
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
            sph.exportState("./output/");
            write_frame(sph.particles, count);
            ++count;
        }

        //Update progress bar
        ++show_progress;
    }

    return 0;
}
