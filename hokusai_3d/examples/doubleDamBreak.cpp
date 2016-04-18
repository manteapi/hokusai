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
    int resolution = 2.5e6; ///particle number per m3
    System<IISPHSolver> sph(resolution);

    Vec3r  fluidBox1(0.5,1.0,0.5);
    Vec3r  fluidOffset1(0,0,0);
    sph.addParticleBox(fluidOffset1, fluidBox1);

    Vec3r  fluidBox2(0.5,1.0,0.5);
    Vec3r  fluidOffset2(1.0,0.0,0.5);
    sph.addParticleBox(fluidOffset2, fluidBox2);

    Vec3r  securityOffset(1.05*sph.getSmoothingRadius());
    Vec3r  boundBox(1.5,1.0,1.0);
    boundBox += securityOffset;
    Vec3r  boundOffset = fluidOffset1;
    boundOffset -= securityOffset;
    sph.addBoundaryBox(boundOffset, boundBox);

    sph.init();

    sph.getViscosity() = 0.01;
    sph.getFluidCohesion() = 0.05;
    sph.getBoundaryAdhesion() = 0.01;
    sph.getBoundaryFriction() = 0.01;
    sph.getTimeStep() = 1e-4;

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
        sph.computeSimulationStep();

        //Output
        if( std::floor((sph.getTime()-sph.getTimeStep())/0.016) != std::floor(sph.getTime()/0.016) )
        {
            write_frame< System<IISPHSolver>::Particle >(sph.getParticles(), count, 10.0);
            sph.exportState("./");
            ++count;
        }

        //Update progress bar
        ++show_progress;
    }

    return 0;
}
