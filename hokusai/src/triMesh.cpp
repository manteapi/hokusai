#ifndef TRIMESH_CPP
#define TRIMESH_CPP

#include "./../include/hokusai/triMesh.hpp"
#include "./../include/hokusai/io.hpp"

TriMesh::TriMesh()
{
    vertices.clear();
    triangles.clear();
    normals.clear();
}

TriMesh::TriMesh(const std::vector< Vec3f >& _vertices, const std::vector< Vec3f >& _normals, const std::vector< Vec3i >& _triangles)
{
    vertices = _vertices;
    normals = _normals;
    triangles = _triangles;
}

TriMesh::TriMesh(const TriMesh& _triMesh)
{
    vertices = _triMesh.vertices;
    normals = _triMesh.normals;
    triangles = _triMesh.triangles;
}

TriMesh::TriMesh(const char* filename)
{
    FILE* pFile=NULL;
    pFile = fopen (filename , "r");
    if (pFile == NULL)
    {
        std::cout << "Error opening file" << std::endl;
    }
    else
    {
        read_obj(pFile, vertices, normals, triangles);
        fclose (pFile);
    }
}

TriMesh::~TriMesh(){}

#endif
