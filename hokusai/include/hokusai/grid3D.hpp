#ifndef GRID3D_H
#define GRID3D_H

#include <array>
#include <vector>
#include <cmath>
#include <iostream>
//#include <GL/gl.h>
#include <algorithm>
#include "Vec.hpp"
#include "particle.hpp"

using namespace std;

class Cell
{
    typedef Vec3<double> Vec;
    public :
        vector< int > fluid;
        vector< int > boundary;

    public :
        Cell()
        {
            fluid.clear(); 
            boundary.clear();
        }

        Cell(const Cell& c)
        {
            fluid = c.fluid;
            boundary = c.boundary;
        }

        ~Cell(){}

    public :
        inline void removeFluid( int i )
        { 
            std::vector<int>::iterator it;
            it = find (fluid.begin(), fluid.end(), i);
            if(it==fluid.end())
            {
                std::cout << i << "Not found" << std::endl; 
                std::cout << *it << std::endl;
                exit(EXIT_FAILURE);
            }
            fluid.erase(it);
        }
        inline void addFluid( int i ){ fluid.push_back(i); }
        inline void addBoundary( int i ){ boundary.push_back(i); }

        void clearFluid(){ fluid.clear(); }
        void clearBoundary(){ boundary.clear(); }
        void clear(){ clearFluid(); clearBoundary(); }

        int fluidNumber(){ return fluid.size(); }
        int boundaryNumber(){ return boundary.size(); }

        vector< int > getFluid(){ return fluid; }
        vector< int > getBoundary(){ return boundary; }
};

class Grid3D
{
    //typedef array<double,3> Vec;
    public :
        typedef Vec3<double> Vec;
        Grid3D();
        Grid3D( const Grid3D& grid );
        ~Grid3D();

        int cellNumber, wNumber, hNumber, dNumber;
        double cellSize, w, h, d;
        Vec minCorner, maxCorner;
        vector< Cell > cells;

    public :
        void init();

        int getCellId( Vec& position );
        array<int,3> getGridCoordinate( Vec& position );
        int addFluid(int id, Particle& p);
        int addFluid(int id, Vec& position);
        int addBoundary(int id, Vec& position);

        vector<int> getNearestFluidNeighbor( Vec & x);
        vector<int> getNearestBoundaryNeighbor( Vec & x);
        void getNearestNeighbor( Vec & x, vector<int>& fluidNeighbor, vector<int>& boundaryNeighbor);
        void getNearestNeighbor( int i, vector<Particle>& particles, vector<Boundary>& boundaries, double radius);
        vector<int> getNearestFluidNeighbor(array<double,3>& x, vector<Particle>& particles, double radius);

        void clearFluid();
        void clearBoundary();
        void clear();

        int fluidNumber();
        int boundaryNumber();

        //void draw();
        void debugCells();
        void info();
};

typedef Vec3<double> Vec;
bool testGridCreation();
bool testSearchNeighbor();
void testAll();

#endif //GRID3D_H
