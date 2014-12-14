#ifndef TRIMESH_HPP
#define TRIMESH_HPP

#include <vector>
#include <array>
#include "Vec.hpp"

typedef Vec3<float> Vec3f;
typedef Vec3<int> Vec3i;

class TriMesh
{
    public :
        TriMesh();
        TriMesh(const std::vector< Vec3f >& _vertices, const std::vector< Vec3f >& _normals, const std::vector< Vec3i>& _triangles);
        TriMesh(const TriMesh& _triMesh);
        TriMesh(const char* filename);
        ~TriMesh();

        std::vector< Vec3f > vertices;
        std::vector< Vec3f > normals;
        std::vector< Vec3i> triangles;
};

#endif
