/*
* Copyright 2015 Pierre-Luc Manteaux
*
*This file is part of Hokusai.
*
*Hokusai is free software: you can redistribute it and/or modify
*it under the terms of the GNU General Public License as published by
*the Free Software Foundation, either version 3 of the License, or
*(at your option) any later version.
*
*Hokusai is distributed in the hope that it will be useful,
*but WITHOUT ANY WARRANTY; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*GNU General Public License for more details.
*
*You should have received a copy of the GNU General Public License
*along with Hokusai.  If not, see <http://www.gnu.org/licenses/>.
*
* Author : Pierre-Luc Manteaux
* Contact : pierre-luc.manteaux@inria.fr
*/

#ifndef HOKUSAI_TRIMESH_CPP
#define HOKUSAI_TRIMESH_CPP

#include "./../include/hokusai/triMesh.hpp"
#include "./../include/hokusai/io.hpp"
#include <set>

namespace hokusai
{

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

void TriMesh::info()
{
    std::cout << "Triangle count : " << triangles.size() << std::endl;
    std::cout << "Vertices count : " << vertices.size() << std::endl;
}

void TriMesh::addCubeMesh(const Vec3f& defaultPosition, const double& scale)
{
    int vertexId;

    //Front face
    vertexId = vertices.size();
    vertices.push_back( defaultPosition + Vec3f(0,0,scale) );
    vertices.push_back( defaultPosition + Vec3f(scale,0,scale) );
    vertices.push_back( defaultPosition + Vec3f(0,scale,scale) );
    triangles.push_back( Vec3i(vertexId, vertexId+1, vertexId+2) );

    vertexId = vertices.size();
    vertices.push_back( defaultPosition + Vec3f(scale,0,scale) );
    vertices.push_back( defaultPosition + Vec3f(0,scale,scale) );
    vertices.push_back( defaultPosition + Vec3f(scale,scale,scale) );
    triangles.push_back( Vec3i(vertexId, vertexId+1, vertexId+2) );

    //Right face
    vertexId = vertices.size();
    vertices.push_back( defaultPosition + Vec3f(scale,scale,scale) );
    vertices.push_back( defaultPosition + Vec3f(scale,0,scale) );
    vertices.push_back( defaultPosition + Vec3f(scale,scale,0) );
    triangles.push_back( Vec3i(vertexId, vertexId+1, vertexId+2) );

    vertexId = vertices.size();
    vertices.push_back( defaultPosition + Vec3f(scale,0,scale) );
    vertices.push_back( defaultPosition + Vec3f(scale,scale,0) );
    vertices.push_back( defaultPosition + Vec3f(scale,0,0) );
    triangles.push_back( Vec3i(vertexId, vertexId+1, vertexId+2) );

    //Back face
    vertexId = vertices.size();
    vertices.push_back( defaultPosition + Vec3f(scale,0,0) );
    vertices.push_back( defaultPosition + Vec3f(0,0,0) );
    vertices.push_back( defaultPosition + Vec3f(scale,scale,0) );
    triangles.push_back( Vec3i(vertexId, vertexId+1, vertexId+2) );

    vertexId = vertices.size();
    vertices.push_back( defaultPosition + Vec3f(0,0,0) );
    vertices.push_back( defaultPosition + Vec3f(scale,scale,0) );
    vertices.push_back( defaultPosition + Vec3f(0,scale,0) );
    triangles.push_back( Vec3i(vertexId, vertexId+1, vertexId+2) );

    //Left face
    vertexId = vertices.size();
    vertices.push_back( defaultPosition + Vec3f(0,scale,0) );
    vertices.push_back( defaultPosition + Vec3f(0,0,0) );
    vertices.push_back( defaultPosition + Vec3f(0,scale,scale) );
    triangles.push_back( Vec3i(vertexId, vertexId+1, vertexId+2) );

    vertexId = vertices.size();
    vertices.push_back( defaultPosition + Vec3f(0,0,0) );
    vertices.push_back( defaultPosition + Vec3f(0,scale,scale) );
    vertices.push_back( defaultPosition + Vec3f(0,0,scale) );
    triangles.push_back( Vec3i(vertexId, vertexId+1, vertexId+2) );

    //Bottom face
    vertexId = vertices.size();
    vertices.push_back( defaultPosition + Vec3f(0,0,scale) );
    vertices.push_back( defaultPosition + Vec3f(0,0,0) );
    vertices.push_back( defaultPosition + Vec3f(scale,0,scale) );
    triangles.push_back( Vec3i(vertexId, vertexId+1, vertexId+2) );

    vertexId = vertices.size();
    vertices.push_back( defaultPosition + Vec3f(0,0,0) );
    vertices.push_back( defaultPosition + Vec3f(scale,0,scale) );
    vertices.push_back( defaultPosition + Vec3f(scale,0,0) );
    triangles.push_back( Vec3i(vertexId, vertexId+1, vertexId+2) );

    //Top face
    vertexId = vertices.size();
    vertices.push_back( defaultPosition + Vec3f(0,scale,scale) );
    vertices.push_back( defaultPosition + Vec3f(scale,scale,scale) );
    vertices.push_back( defaultPosition + Vec3f(0,scale,0) );
    triangles.push_back( Vec3i(vertexId, vertexId+1, vertexId+2) );

    vertexId = vertices.size();
    vertices.push_back( defaultPosition + Vec3f(scale,scale,scale) );
    vertices.push_back( defaultPosition + Vec3f(0,scale,0) );
    vertices.push_back( defaultPosition + Vec3f(scale,scale,0) );
    triangles.push_back( Vec3i(vertexId, vertexId+1, vertexId+2) );
}

void TriMesh::getVertexNeighbor(const int vertexId, std::vector<int>& neighbors)
{
    neighbors.clear();
    for(size_t i=0; i<triangles.size(); ++i)
    {
        for(int j=0; j<3; ++j)
        {
            if(triangles[i][j]==vertexId)
            {
                neighbors.push_back( (triangles[i][(j+1)%3]));
                neighbors.push_back( (triangles[i][(j+2)%3]));
            }
        }
    }
}

void TriMesh::getVertexAdjacentFaces(const int vertexId, std::vector<int>& neighbors)
{
    neighbors.clear();
    std::set<int> tmpNeighborFaces;
    for(size_t i=0; i<triangles.size(); ++i)
    {
        for(int j=0; j<3; ++j)
        {
            if(triangles[i][j]==vertexId)
            {
                tmpNeighborFaces.insert(i);
            }
        }
    }
    for (std::set<int>::iterator it=tmpNeighborFaces.begin(); it!=tmpNeighborFaces.end(); ++it)
    {
        neighbors.push_back(*it);
    }
}

void TriMesh::getAcrossFaces(const int faceId, std::vector<int>& neighbors)
{
    std::set<int> tmpNeighborFaces;
    Vec3i vertex = triangles[faceId];
    std::vector< int > neighborFaces;

    for(int l=0; l<3; ++l)
    {
        getVertexAdjacentFaces(vertex[l], neighborFaces);
        for(size_t i=0; i<neighborFaces.size(); ++i)
        {
            int commonVertex = 0;
            for(int j=0; j< 3; ++j)
            {
                for(int k=0; k<3; ++k)
                {
                    if(triangles[neighborFaces[i]][j]==vertex[k])
                    {
                        commonVertex++;
                    }
                }
            }
            if(commonVertex==2)
            {
                tmpNeighborFaces.insert(neighborFaces[i]);
            }
        }
    }

    neighbors.clear();
    for (std::set<int>::iterator it=tmpNeighborFaces.begin(); it!=tmpNeighborFaces.end(); ++it)
    {
        neighbors.push_back(*it);
    }
}

void TriMesh::getEdges(const std::vector<Vec3i>& _triangles, std::vector< std::pair<int, int> >& _edges) const
{
    std::set< std::pair<int,int> > tmpCleaningBag;
    std::pair<std::set< std::pair<int,int> >::iterator,bool> ret1, ret2;
    for(size_t i=0; i<_triangles.size(); ++i)
    {
        ret1 = tmpCleaningBag.insert( std::pair<int, int>( _triangles[i][0], _triangles[i][1] ) );
        ret2 = tmpCleaningBag.insert( std::pair<int, int>( _triangles[i][1], _triangles[i][0] ) );
        if(ret1.second == true && ret2.second == true)
            tmpCleaningBag.erase(ret2.first);
        if (ret1.second == false && ret2.second == true)
            tmpCleaningBag.erase(ret2.first);
        if (ret1.second == true && ret2.second == false)
            tmpCleaningBag.erase(ret1.first);


        ret1 = tmpCleaningBag.insert( std::pair<int, int>( _triangles[i][1], _triangles[i][2] ) );
        ret2 = tmpCleaningBag.insert( std::pair<int, int>( _triangles[i][2], _triangles[i][1] ) );
        if(ret1.second == true && ret2.second == true)
            tmpCleaningBag.erase(ret2.first);
        if (ret1.second == false && ret2.second == true)
            tmpCleaningBag.erase(ret2.first);
        if (ret1.second == true && ret2.second == false)
            tmpCleaningBag.erase(ret1.first);

        ret1 = tmpCleaningBag.insert( std::pair<int, int>( _triangles[i][2], _triangles[i][0] ) );
        ret2 = tmpCleaningBag.insert( std::pair<int, int>( _triangles[i][0], _triangles[i][2] ) );
        if(ret1.second == true && ret2.second == true)
            tmpCleaningBag.erase(ret2.first);
        if (ret1.second == false && ret2.second == true)
            tmpCleaningBag.erase(ret2.first);
        if (ret1.second == true && ret2.second == false)
            tmpCleaningBag.erase(ret1.first);
    }
    for(std::set< std::pair<int,int> >::iterator it=tmpCleaningBag.begin(); it!=tmpCleaningBag.end(); ++it)
        _edges.push_back(*it);
}

void TriMesh::getBorderEdges(const int triangleId, std::vector< std::pair<int,int> >& edges)
{
    edges.clear();
    std::pair<int,int> tmpEdge;
    std::vector<int> tmpTriangles;
    for(int i=0; i<3; ++i)
    {
        tmpEdge.first = triangles[triangleId][i%3];
        tmpEdge.second = triangles[triangleId][(i+1)%3];
        tmpTriangles = adjacentFaces[tmpEdge.first];
        int acrossEdgeTriangle = 0;
        for(size_t j=0; j<tmpTriangles.size(); ++j)
        {
            for(int k=0; k<3; ++k)
            {
                if(triangles[tmpTriangles[j]][k]==tmpEdge.second)
                    acrossEdgeTriangle++;
            }
        }
        if(acrossEdgeTriangle==1)
            edges.push_back(tmpEdge);
    }
}

//To do : debug
void TriMesh::getContour(std::vector< std::pair<int, int> >& /*edges*/)
{
    /*
    edges.clear();

    //Look for one border triangle
    int borderTriangle=-1;
    for(size_t i=0; i<triangles.size(); ++i)
    {
        if(acrossFaces[i].size() != 3)
        {
            borderTriangle = i;
            break;
        }
    }
    if(borderTriangle == -1){ return; }

    //Look for one border edge of the border triangle
    std::vector< std::pair<int,int> > tmpEdges;

    //Follow the border edge from border to border until going back home
    int tmp;
    std::pair<int,int> startEdge, currentEdge, nextEdge;
    std::vector<int> tmpTriangles;
    startEdge = tmpEdges[0];
    currentEdge = startEdge;
    edges.push_back(startEdge);
    while( nextEdge != startEdge )
    {
        bool foundNextEdge = false;
        tmpTriangles = adjacentFaces[currentEdge.second];
        for(int i=0; i<tmpTriangles.size(); ++i)
        {
            getBorderEdges(tmpTriangles[i], tmpEdges);
            for(int j=0; j<tmpEdges.size(); ++j)
            {
                if(tmpEdges[j].first==currentEdge.first || tmpEdges[j].second==currentEdge.second)
                {
                    if(tmpEdges[j].first!=currentEdge.first || tmpEdges[j].second!=currentEdge.second)
                    {
                        nextEdge = tmpEdges[j];
                        foundNextEdge;
                        break;
                    }
                }
            }
            if(foundNextEdge == true) break;
        }

        //Put the next edge in the right direction
        if(nextEdge.first != currentEdge.second)
        {
            tmp = nextEdge.first;
            nextEdge.first = nextEdge.second;
            nextEdge.second = tmp;
        }

        edges.push_back(nextEdge);
        currentEdge = nextEdge;
    }
    */
}

void TriMesh::clean()
{
    vertices.clear();
    triangles.clear();
    normals.clear();
    vneighbors.clear();
    adjacentFaces.clear();
    acrossFaces.clear();
}

void TriMesh::buildNeighbors()
{

    vneighbors.clear();
    adjacentFaces.clear();
    vneighbors.resize(vertices.size());
    adjacentFaces.resize(vertices.size());

    for(size_t i=0; i<vertices.size(); ++i)
    {
        std::vector<int> tmpVertexNeighbor;
        std::set<int> tmpCleaningBag;
        getVertexNeighbor(i, tmpVertexNeighbor);
        for(size_t j=0; j<tmpVertexNeighbor.size(); ++j)
            tmpCleaningBag.insert(tmpVertexNeighbor[j]);
        for(std::set<int>::iterator it=tmpCleaningBag.begin(); it!=tmpCleaningBag.end(); ++it)
            vneighbors[i].push_back(*it);

        std::vector<int> tmpAdjacentFaces;
        getVertexAdjacentFaces(i, tmpAdjacentFaces);
        adjacentFaces[i] = tmpAdjacentFaces;
    }


    acrossFaces.clear();
    acrossFaces.resize(triangles.size());
    for(size_t i=0; i<triangles.size(); ++i)
    {
        std::vector<int> tmpAcrossFaces;
        getAcrossFaces(i, tmpAcrossFaces);
        acrossFaces[i] = tmpAcrossFaces;
    }
}
}
#endif
