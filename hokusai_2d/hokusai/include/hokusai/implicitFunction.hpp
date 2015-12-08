#ifndef IMPLICIT_FUNCTION_HPP
#define IMPLICIT_FUNCTION_HPP

#include "utility.hpp"
#include "gridUtility.hpp"
#include "system.hpp"
#include <limits.h>

namespace hokusai
{

void computeScalarField(std::vector<double>& scalarField, Grid2dUtility& gridInfo, System &fluid,
                        const double& resolution, const double &initialValue, const double& radius);
double implicit_sphere(const Vec2d & x1, const Vec2d & x2, const double radius);
double implicit_metaball(const Vec2d & x1, const Vec2d & x2, const double radius);

} // namespace hokusai

#endif //IMPLICIT_FUNCTION_HPP
