#ifndef SAMPLER_HPP
#define SAMPLER_HPP

#include <aljabr/Vec.hpp>
#include <shkodra/tri2mesh.hpp>
#include <vector>

typedef double Real;
typedef aljabr::Vec3<Real> Vec3r;
typedef aljabr::Vec2<Real> Vec2r;
typedef shkodra::Tri2Mesh Mesh;

bool LineLineIntersect(const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const Vec3r& p4, Vec3r& pa, Vec3r& pb,float & mua, float & mub);
bool LineIntersect(const Vec2r& p1, const Vec2r& p2, const Vec2r& p3, const Vec2r& p4, Vec2r& p);
bool LineIntersect(const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const Vec3r& p4, Vec3r& p);
bool AkinciFullTriangleSampling( const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const float& particleDiameter, std::vector< Vec3r >& samples);
bool AkinciTriangleSampling( const Vec3r& p1, const Vec3r& p2, const Vec3r& p3, const float& particleDiameter, std::vector< Vec3r >& samples);
bool AkinciEdgeSampling( const Vec3r& p1, const Vec3r& p2, const float& particleDiameter, std::vector< Vec3r >& samples);
bool AkinciMeshSampling(const Mesh& mesh, const float &particleDiameter, std::vector<Vec3r>& samples);
#endif
