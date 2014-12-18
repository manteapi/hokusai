#ifndef HOKUSAI_SAMPLER_HPP
#define HOKUSAI_SAMPLER_HPP

#include "Vec.hpp"
#include "triMesh.hpp"
#include <vector>

namespace hokusai
{
bool LineLineIntersect(const Vec3<double>& p1, const Vec3<double>& p2, const Vec3<double>& p3, const Vec3<double>& p4, Vec3<double>& pa, Vec3<double>& pb,double & mua, double & mub);
bool LineIntersect(const Vec2<double>& p1, const Vec2<double>& p2, const Vec2<double>& p3, const Vec2<double>& p4, Vec2<double>& p);
bool LineIntersect(const Vec3<double>& p1, const Vec3<double>& p2, const Vec3<double>& p3, const Vec3<double>& p4, Vec3<double>& p);
bool AkinciTriangleSampling( const Vec3<double>& p1, const Vec3<double>& p2, const Vec3<double>& p3, const double& particleDiameter, std::vector< Vec3<double> >& samples);
bool AkinciEdgeSampling( const Vec3<double>& p1, const Vec3<double>& p2, const double& particleDiameter, std::vector< Vec3<double> >& samples);
bool AkinciMeshSampling(const TriMesh& mesh, const double &particleDiameter, std::vector<Vec3f>& samples);
}
#endif
