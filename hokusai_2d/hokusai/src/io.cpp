#include "./../include/hokusai/io.hpp"

#include <iostream>
#include <fstream>

namespace hokusai
{

void exportOBJ(const std::string & filename, const std::vector<Edge> & edges)
{
    std::ofstream fileStream;
    fileStream.open(filename.c_str());
    for(size_t i=0; i<edges.size(); ++i)
    {
        fileStream << "v " << edges[i].p1[0] << " " << edges[i].p1[1] << " 0.0" << std::endl;
        fileStream << "v " << edges[i].p2[0] << " " << edges[i].p2[1] << " 0.0" << std::endl;
    }
    for(size_t i=0; i<edges.size(); ++i)
    {
        fileStream << "l " << 2*i+1 << " " << 2*i+2 << std::endl;
    }
    fileStream.close();
}

}//namespace hokusai
