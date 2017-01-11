#include <hokusai/system.hpp>
#include <hokusai/io.hpp>

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
    std::string filename1 = "test_000000_00.bphys";
    std::string filename2 = "test_000001_00.bphys";
    BlenderExporter debugger;
    debugger.debugHeader(filename1);
    debugger.debugFrame(filename2);
}
