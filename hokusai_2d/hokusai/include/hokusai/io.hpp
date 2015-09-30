#ifndef HOKUSAI_IO_H
#define HOKUSAI_IO_H

#include <vector>
#include <string>
#include "marchingSquare.hpp"

namespace hokusai
{

void exportOBJ(const std::string & filename, const std::vector<Edge> & edges);

}

#endif
