#ifndef HOKUSAI_TRIMESH_HPP
#define HOKUSAI_TRIMESH_HPP

#include <vector>
#include <array>
#include <set>
#include <utility>
#include "common.hpp"

namespace hokusai
{
    class TriMesh
    {
        public :
            TriMesh();
            TriMesh(const std::vector< Vec3r >& _vertices, const std::vector< Vec3r >& _normals, const std::vector< Vec3i>& _triangles);
            TriMesh(const TriMesh& _triMesh);
            TriMesh(const char* filename);
            ~TriMesh();

            std::vector< Vec3r > vertices;
            std::vector< Vec3r > normals;
            std::vector< Vec3i> triangles;

            std::vector< std::vector<int> > vneighbors;
            std::vector< std::vector<int> > adjacentFaces;
            std::vector< std::vector<int> > acrossFaces;

            void addCubeMesh(const Vec3r& defaultPosition, const HReal& scale);

            void getVertexNeighbor(const int vertexId, std::vector<int>& vneighbors);
            void getVertexAdjacentFaces(const int vertexId, std::vector<int>& vneighbors);
            void getAcrossFaces(const int faceId, std::vector<int>& vneighbors);
            void getEdges(const std::vector<Vec3i>& _triangles, std::vector< std::pair<int,int> >& _edges) const;
            void getBorderEdges(const int triangleId, std::vector< std::pair<int, int> >& edges);

            void getContour(std::vector< std::pair<int,int> >& edges);

            void buildNeighbors();
            void clean();
            void info();
            bool testBuildNeighbors()
            {
                clean();

                vertices.push_back( Vec3r(-1,0,0) );
                vertices.push_back( Vec3r(1,0,0) );
                vertices.push_back( Vec3r(1,0,0) );
                vertices.push_back( Vec3r(1,1,0) );
                vertices.push_back( Vec3r(-1,1,0) );
                vertices.push_back( Vec3r(0,-1,0) );
                vertices.push_back( Vec3r(1,-1,0) );

                triangles.push_back( Vec3i(0,1,2) );
                triangles.push_back( Vec3i(1,2,3) );
                triangles.push_back( Vec3i(0,2,4) );
                triangles.push_back( Vec3i(0,1,5) );
                triangles.push_back( Vec3i(1,5,6) );

                buildNeighbors();

                bool success = true;
                //Check vneighbors
                success = success && (vertices.size() == 7);
                success = success && (triangles.size() == 5);
                if(!success) return success;


                //Check number of vneighbors
                success = success && ( vneighbors[0].size()==4 ) && ( vneighbors[1].size()==5 )
                    && ( vneighbors[2].size() == 4) && (vneighbors[3].size() == 2)
                    && (vneighbors[4].size() == 2) && (vneighbors[5].size() == 3)
                    && (vneighbors[6].size() == 2);
                if(!success) return success;


                //Check vertex neighbors
                success = success && (vneighbors[0][0] == 1) && (vneighbors[0][1] == 2) && (vneighbors[0][2] == 4) && (vneighbors[0][3] == 5);
                success = success && (vneighbors[1][0] == 0) && (vneighbors[1][1] == 2) && (vneighbors[1][2] == 3) && (vneighbors[1][3] == 5) && (vneighbors[1][4] == 6);
                success = success && (vneighbors[2][0] == 0) && (vneighbors[2][1] == 1) && (vneighbors[2][2] == 3) && (vneighbors[2][3] == 4);
                success = success && (vneighbors[3][0] == 1) && (vneighbors[3][1] == 2);
                success = success && (vneighbors[4][0] == 0) && (vneighbors[4][1] == 2);
                success = success && (vneighbors[5][0] == 0) && (vneighbors[5][1] == 1);

                //Check adjacency faces
                success = success && (adjacentFaces.size() == 7);
                if(!success) return false;

                success = success && ( adjacentFaces[0].size()==3 ) && ( adjacentFaces[1].size()==4 )
                    && ( adjacentFaces[2].size() == 3) && (adjacentFaces[3].size() == 1)
                    && (adjacentFaces[4].size() == 1) && (adjacentFaces[5].size() == 2)
                    && (adjacentFaces[6].size() == 1);
                if(!success) return false;



                success = success && (adjacentFaces[0][0] == 0) && (adjacentFaces[0][1] == 2) && (adjacentFaces[0][2] == 3);
                success = success && (adjacentFaces[1][0] == 0) && (adjacentFaces[1][1] == 1) && (adjacentFaces[1][2] == 3) && (adjacentFaces[1][3] == 4);
                success = success && (adjacentFaces[2][0] == 0) && (adjacentFaces[2][1] == 1) && (adjacentFaces[2][2] == 2);
                success = success && (adjacentFaces[3][0] == 1);
                success = success && (adjacentFaces[4][0] == 2);
                success = success && (adjacentFaces[5][0] == 3);

                //Check accross faces
                success = success && (acrossFaces.size() == 5);
                if(!success) return false;

                success = success && ( acrossFaces[0].size()==3 ) && ( acrossFaces[1].size()==1 )
                    && ( acrossFaces[2].size() == 1) && (acrossFaces[3].size() == 2) && (acrossFaces[4].size() == 1);
                if(!success) return false;

                success = success && (acrossFaces[0][0] == 1) && (acrossFaces[0][1] == 2) && (acrossFaces[0][2] == 3);
                success = success && (acrossFaces[1][0] == 0);
                success = success && (acrossFaces[2][0] == 0);
                success = success && (acrossFaces[3][0] == 0) && (acrossFaces[3][1] == 4);
                success = success && (acrossFaces[4][0] == 3);

                return success;
            }

            bool testGetEdges()
            {
                clean();

                bool success = true;

                triangles.push_back( Vec3i(0,1,2) );
                triangles.push_back( Vec3i(1,2,3) );
                triangles.push_back( Vec3i(0,2,4) );

                std::vector< std::pair<int,int> > edges;
                getEdges(triangles, edges);

                success = success && (edges.size()==7);
                success = success && (edges[0].first == 0 && edges[0].second == 1);
                success = success && (edges[1].first == 1 && edges[1].second == 2);
                success = success && (edges[2].first == 2 && edges[2].second == 0);
                success = success && (edges[3].first == 2 && edges[3].second == 3);
                success = success && (edges[4].first == 2 && edges[4].second == 4);
                success = success && (edges[5].first == 3 && edges[5].second == 1);
                success = success && (edges[6].first == 4 && edges[6].second == 0);

                return success;

            }

            bool testGetBorderEdges()
            {
                clean();

                bool success = true;

                vertices.push_back( Vec3r(-1,0,0) );
                vertices.push_back( Vec3r(1,0,0) );
                vertices.push_back( Vec3r(1,0,0) );
                vertices.push_back( Vec3r(1,1,0) );
                vertices.push_back( Vec3r(-1,1,0) );
                vertices.push_back( Vec3r(0,-1,0) );
                vertices.push_back( Vec3r(1,-1,0) );

                triangles.push_back( Vec3i(0,1,2) );
                triangles.push_back( Vec3i(1,2,3) );
                triangles.push_back( Vec3i(0,2,4) );
                triangles.push_back( Vec3i(0,1,5) );
                triangles.push_back( Vec3i(1,5,6) );

                buildNeighbors();

                std::vector< std::pair<int,int> > borderEdge;

                getBorderEdges(0, borderEdge);
                success = success && (borderEdge.size() == 0);
                getBorderEdges(1, borderEdge);
                success = success && (borderEdge.size() == 2);
                getBorderEdges(2, borderEdge);
                success = success && (borderEdge.size() == 2);
                getBorderEdges(3, borderEdge);
                success = success && (borderEdge.size() == 1);
                getBorderEdges(4, borderEdge);
                success = success && (borderEdge.size() == 2);

                getBorderEdges(1, borderEdge);
                success = success && (borderEdge[0].first== 2) && (borderEdge[0].second==3);
                success = success && (borderEdge[1].first== 3) && (borderEdge[1].second==1);

                getBorderEdges(2, borderEdge);
                success = success && (borderEdge[0].first== 2) && (borderEdge[0].second==4);
                success = success && (borderEdge[1].first== 4) && (borderEdge[1].second==0);

                getBorderEdges(3, borderEdge);
                success = success && (borderEdge[0].first== 5) && (borderEdge[0].second==0);

                getBorderEdges(4, borderEdge);
                success = success && (borderEdge[0].first== 5) && (borderEdge[0].second==6);
                success = success && (borderEdge[1].first== 6) && (borderEdge[1].second==1);

                return success;
            }

            bool testGetContour()
            {
                bool success = true;

                vertices.push_back( Vec3r(-1,0,0) );
                vertices.push_back( Vec3r(1,0,0) );
                vertices.push_back( Vec3r(1,0,0) );
                vertices.push_back( Vec3r(1,1,0) );
                vertices.push_back( Vec3r(-1,1,0) );
                vertices.push_back( Vec3r(0,-1,0) );
                vertices.push_back( Vec3r(1,-1,0) );

                triangles.push_back( Vec3i(0,1,2) );
                triangles.push_back( Vec3i(1,2,3) );
                triangles.push_back( Vec3i(0,2,4) );
                triangles.push_back( Vec3i(0,1,5) );
                triangles.push_back( Vec3i(1,5,6) );

                buildNeighbors();

                std::vector< std::pair<int, int> > edges;
                getContour(edges);
                return success;
            }
    };
}
#endif
