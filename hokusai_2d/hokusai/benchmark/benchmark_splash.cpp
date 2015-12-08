#include <iostream>
#include <hokusai/system.hpp>

using namespace std;

int main()
{
    int counter = 0, maxIter = 1e3;
    int resolution = 4e3;
    double timeStep = 2e-3;
    hokusai::System sph(resolution);
    sph.set_benchmark_splash(timeStep);
    sph.init();
    while(counter<maxIter)
    {
        counter++;
        sph.simulate();
    }
    std::cout << "Benchmark dam break done" << std::endl;

    return 0;
}
