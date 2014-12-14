#include "../include/hokusai/grid3D.hpp"

using namespace std;

Grid3D::Grid3D()
{
    cellNumber = wNumber = hNumber = dNumber = 0;
    cellSize = w = h = d = 0.0;
    minCorner = maxCorner = Vec(0.0f);
    cells.clear();
}

Grid3D::Grid3D(const Grid3D& grid)
{
    cellNumber = grid.cellNumber;
    wNumber = grid.wNumber;
    hNumber = grid.hNumber;
    dNumber = grid.dNumber;

    cellSize = grid.cellSize;
    w = grid.w;
    h = grid.h;
    d = grid.d;
    minCorner = grid.minCorner;
    maxCorner = grid.maxCorner;
    cells = grid.cells;
}

Grid3D::~Grid3D(){}

void Grid3D::init()
{
    //First shot
    w = maxCorner[0] - minCorner[0];
    h = maxCorner[1] - minCorner[1];
    d = maxCorner[2] - minCorner[2];
    if(cellSize <= 0){std::cout << "cellsize uninitialized" << std::endl;}
    wNumber = w/(double)cellSize;
    hNumber = h/(double)cellSize;
    dNumber = d/(double)cellSize;

    //Correction
    maxCorner[0] = minCorner[0] + wNumber*cellSize;
    maxCorner[1] = minCorner[1] + hNumber*cellSize;
    maxCorner[2] = minCorner[2] + dNumber*cellSize;
    w = maxCorner[0] - minCorner[0];
    h = maxCorner[1] - minCorner[1];
    d = maxCorner[2] - minCorner[2];
    wNumber = w/(double)cellSize;
    hNumber = h/(double)cellSize;
    dNumber = d/(double)cellSize;
    cellNumber = wNumber*hNumber*dNumber;
    cells.resize(cellNumber);

}

int Grid3D::addFluid(int id, Particle& p)
{
    array<int,3> coord = getGridCoordinate(p.x);
    int cellId = coord[0] + coord[1]*wNumber+ coord[2]*wNumber*hNumber;

    if( cellId < 0 || cellId >= (int)cells.size() )
    {
        std::cerr << "addFluid : cell " << cellId << " out of range" << std::endl;
        return -1;
    }
    cells[cellId].addFluid(id);
    return cellId;
}

int Grid3D::addFluid(int id, Vec& position)
{
    array<int,3> coord = getGridCoordinate(position);
    //int cellId = coord[2] + coord[1]*hNumber+ coord[0]*wNumber*hNumber;
    int cellId = coord[0] + coord[1]*wNumber+ coord[2]*wNumber*hNumber;
    if( cellId < 0 || cellId >= (int)cells.size() )
    {
        std::cerr << "addFluid : cell " << cellId << " out of range" << std::endl;
        return -1;
    }
    cells[cellId].addFluid(id);
    return cellId;
}

int Grid3D::getCellId( Vec& position )
{
    array<int,3> coord = getGridCoordinate(position);
    int cellId = coord[0] + coord[1]*wNumber+ coord[2]*wNumber*hNumber;
    return cellId;
}

array<int,3> Grid3D::getGridCoordinate( Vec& position )
{
    int i = std::floor( (position[0]-minCorner[0])/(double)cellSize ); //width
    int j = std::floor( (position[1]-minCorner[1])/(double)cellSize ); //height
    int k = std::floor( (position[2]-minCorner[2])/(double)cellSize ); //depth
    array<int,3> index = {{i,j,k}};
    return index;
}

int Grid3D::addBoundary(int id, Vec& position)
{
    array<int,3> coord = getGridCoordinate(position);
    int cellId = coord[0] + coord[1]*wNumber+ coord[2]*wNumber*hNumber;
    if( cellId < 0 || cellId >= (int)cells.size() )
    {
        std::cerr << "addBoundary : cell "<< cellId << " : out of range" << std::endl;
        std::cout << "Position : " << position[0] << ", " << position[1] << " , " << position[2] << std::endl;
        std::cout << "Coord : " << coord[0] << ", " << coord[1] << " , " << coord[2] << std::endl;
        std::cout << "w,h,d : " << w << ", " << h << ", " << d << std::endl;
        std::cout << "number : " << cells.size() << std::endl;
        std::cout << "cell number : " << cellNumber << std::endl;
        std::cout << "Max : " << maxCorner[0] << ", " << maxCorner[1] << ", " << maxCorner[2] << std::endl;
        std::cout << "Min : " << minCorner[0] << ", " << minCorner[1] << ", " << minCorner[2] << std::endl;
        return -1;
    }
    cells[cellId].addBoundary(id);
    return cellId;
}

int Grid3D::fluidNumber()
{
    int number = 0;
    for(Cell& c_i : cells)
        number += c_i.fluidNumber();
    return number;
}

int Grid3D::boundaryNumber()
{
    int number = 0;
    for(Cell& c_i : cells)
        number += c_i.boundaryNumber();
    return number;
}

void Grid3D::info()
{
    std::cout << "Cell Number : " << cellNumber << std::endl;
    std::cout << "Cell Number : " << wNumber << " x " << hNumber << " x " << dNumber << std::endl;
    std::cout << "Dimensions : " << w << " x " << h << " x " << d << std::endl;
    std::cout << "Cell Size : " << cellSize << std::endl;
}

void Grid3D::debugCells()
{
    int i = 0;
    for(Cell& c : cells)
    {
        std::cout << "Cell : " << i << std::endl;
        ++i;
        for(int& id : c.fluid)
        {
            std::cout << "Fluid id : " << id << std::endl;
        }
    }
}

/*
void Grid3D::draw()
{
    //Get grid size
    double width = maxCorner[0] - minCorner[0];
    double height = maxCorner[1] - minCorner[1];
    double depth = maxCorner[2] - minCorner[2];

    int widthSize = std::floor(width/(double)cellSize); //width cell number
    int heightSize = std::floor(height/(double)cellSize); // height cell number
    int depthSize = std::floor(depth/(double)cellSize); // depth cell number

    //Depth Line
    for(int k = 0; k <= depthSize; ++k)
    {
        //Vertical Line
        for(int i =0; i <= widthSize; ++i)
        {
            //glLineWidth(1.0);
            glBegin(GL_LINES);
            glColor3f(1.0,1.0,1.0);
            glVertex3f(minCorner[0]+cellSize*i, minCorner[1], minCorner[2]+cellSize*k);
            glVertex3f(minCorner[0]+cellSize*i, maxCorner[1], minCorner[2]+cellSize*k);
            glEnd();
        }
        //Horizontal Line
        for(int i =0; i <= heightSize; ++i)
        {
            //glLineWidth(0.5);
            glBegin(GL_LINES);
            //glColor3f(0.5,0.5,0.5);
            glColor3f(1.0,1.0,1.0);
            glVertex3f(minCorner[0], minCorner[1]+cellSize*i, minCorner[2]+cellSize*k);
            glVertex3f(maxCorner[0], minCorner[1]+cellSize*i, minCorner[2]+cellSize*k);
            glEnd();
        }
    }

    for(int i = 0; i <= widthSize; ++i)
    {
        for(int j = 0; j <= heightSize; ++j)
        {
            glBegin(GL_LINES);
            //glColor3f(0.5,0.5,0.5);
            glColor3f(1.0,1.0,1.0);
            glVertex3f(minCorner[0]+cellSize*i, minCorner[1]+cellSize*j, minCorner[2]);
            glVertex3f(minCorner[0]+cellSize*i, minCorner[1]+cellSize*j, maxCorner[2]);
            glEnd();
        }
    }
}
*/

void Grid3D::clearFluid()
{
    for(Cell& c : cells)
        c.clearFluid();
}

void Grid3D::clearBoundary()
{
    for(Cell& c : cells)
        c.clearBoundary();
}

void Grid3D::clear()
{
    for(Cell& c : cells)
        c.clear();
}

//vector<int> Grid3D::getNearestFluidNeighbor( Vec & x )
//{
//    //Get grid coordinate
//    array<int,3> coord = getGridCoordinate(x);

//    //Get bounding box coordinate
//    int minI = max(coord[0]-1,0);
//    int maxI = min(coord[0]+1,wNumber-1);

//    int minJ = max(coord[1]-1,0);
//    int maxJ = min(coord[1]+1,hNumber-1);

//    int minK = max(coord[2]-1,0);
//    int maxK = min(coord[2]+1,dNumber-1);

//    vector< int > neighborsId;
//    for( int _d = minK; _d <= maxK; _d++ )
//    {
//        for( int _h = minJ; _h <= maxJ; _h++)
//        {
//            for( int _w = minI; _w <= maxI; _w++)
//            {
//                //Get the cell ID
//                int cellId = _w + _h*wNumber+ _d*wNumber*hNumber;
//                if( cellId < 0 || cellId >= (int)cells.size() )
//                {
//                    std::cerr << "Warning getNearestFluidNeighbor : " << cellId << std::endl;
//                    return neighborsId;
//                }
//                //std::cout << cellId << std::endl;
//                Cell& c = cells[cellId];
//                int size = c.fluidNumber();
//                //std::cout << "Taille : " << size << std::endl;
//                for(int l = 0; l < size; ++l)
//                {
//                    neighborsId.push_back(c.fluid[l]);
//                }
//            }
//        }
//    }
//    return neighborsId;
//}

//vector<int> Grid3D::getNearestBoundaryNeighbor( Vec & x )
//{
//    array<int,3> coord = getGridCoordinate(x);

//    //Get bounding box coordinate
//    vector< int > neighborsId;
//    for( int _d = max(coord[2]-1,0); _d <= min(coord[2]+1,dNumber-1); ++_d )
//    {
//        for( int _h = max(coord[1]-1,0); _h <= min(coord[1]+1,hNumber-1); ++_h)
//        {
//            for( int _w = max(coord[0]-1,0); _w <= min(coord[0]+1,wNumber-1); ++_w)
//            {
//                //Get the cell ID
//                int cellId = _w + _h*wNumber+ _d*wNumber*hNumber;
//                if( cellId < 0 || cellId >= (int)cells.size() )
//                {
//                    std::cerr << "Warning : " << cellId << std::endl;
//                    return neighborsId;
//                }
//                Cell& c = cells[cellId];
//                int size = c.boundaryNumber();
//                for(int l = 0; l < size; ++l)
//                {
//                    neighborsId.push_back(c.boundary[l]);
//                }
//            }
//        }
//    }
//    return neighborsId;
//}

vector<int> Grid3D::getNearestBoundaryNeighbor( Vec& x )
{
    vector<int> boundaryNeighbor;
    array<int,3> coord = getGridCoordinate(x);
    for( int _d = max(coord[2]-1,0); _d <= min(coord[2]+1,dNumber-1); ++_d )
        for( int _h = max(coord[1]-1,0); _h <= min(coord[1]+1,hNumber-1); ++_h)
            for( int _w = max(coord[0]-1,0); _w <= min(coord[0]+1,wNumber-1); ++_w){
                int cellId = _w + _h*wNumber+ _d*wNumber*hNumber;
                if( cellId < 0 || cellId >= (int)cells.size() )
                {
                    std::cerr << "Warning : " << cellId << std::endl;
                    boundaryNeighbor.clear();
                    return boundaryNeighbor;
                }
                Cell& c = cells[cellId];
                for(int& j : c.boundary)
                    boundaryNeighbor.push_back(j);
            }
    return boundaryNeighbor;
}

vector<int> Grid3D::getNearestFluidNeighbor( Vec& x )
{
    vector<int> fluidNeighbor;
    array<int,3> coord = getGridCoordinate(x);
    for( int _d = max(coord[2]-1,0); _d <= min(coord[2]+1,dNumber-1); ++_d )
        for( int _h = max(coord[1]-1,0); _h <= min(coord[1]+1,hNumber-1); ++_h)
            for( int _w = max(coord[0]-1,0); _w <= min(coord[0]+1,wNumber-1); ++_w){
                int cellId = _w + _h*wNumber+ _d*wNumber*hNumber;
                if( cellId < 0 || cellId >= (int)cells.size() )
                {
                    std::cerr << "Warning : " << cellId << std::endl;
                    fluidNeighbor.clear();
                    return fluidNeighbor;
                }
                Cell& c = cells[cellId];
                for(int& j : c.fluid)
                    fluidNeighbor.push_back(j);
            }
    return fluidNeighbor;
}

void Grid3D::getNearestNeighbor( Vec & x, vector<int>& fluidNeighbor, vector<int>& boundaryNeighbor)
{
    array<int,3> coord = getGridCoordinate(x);
    for( int _d = max(coord[2]-1,0); _d <= min(coord[2]+1,dNumber-1); _d++ ){
        for( int _h = max(coord[1]-1,0); _h <= min(coord[1]+1,hNumber-1); _h++){
            for( int _w = max(coord[0]-1,0); _w <= min(coord[0]+1,wNumber-1); _w++){
                int cellId = _w + _h*wNumber+ _d*wNumber*hNumber;
                if( cellId < 0 || cellId >= (int)cells.size() )
                {
                    std::cerr << "Warning : " << cellId << std::endl;
                    return ;
                }
                Cell& c = cells[cellId];
                for(int& j : c.boundary)
                    boundaryNeighbor.push_back(j);
                for(int& j : c.fluid)
                    fluidNeighbor.push_back(j);
            }}}
}

void Grid3D::getNearestNeighbor( int i, vector<Particle>& particles, vector<Boundary>& boundaries, double radiusSquared)
{
    Particle& pi = particles[i];
    pi.fluidNeighbor.clear();
    pi.boundaryNeighbor.clear();
    array<int,3> coord = getGridCoordinate(pi.x);
    for( int _d = max(coord[2]-1,0); _d <= min(coord[2]+1,dNumber-1); ++_d )
        for( int _h = max(coord[1]-1,0); _h <= min(coord[1]+1,hNumber-1); ++_h)
            for( int _w = max(coord[0]-1,0); _w <= min(coord[0]+1,wNumber-1); ++_w){
                int cellId = _w + _h*wNumber+ _d*wNumber*hNumber;
                if( cellId < 0 || cellId >= (int)cells.size() )
                {
                    std::cerr << "Warning : " << cellId << std::endl;
                    return ;
                }

                Cell& c = cells[cellId];
                for(int& j : c.boundary)
                {
                    double dist = (pi.x-boundaries[j].x).lengthSquared();
                    if(dist<radiusSquared)
                        pi.boundaryNeighbor.push_back(j);
                }
                for(int& j : c.fluid)
                {
                    double dist = (pi.x-particles[j].x).lengthSquared();
                    if(dist<=radiusSquared)
                        pi.fluidNeighbor.push_back(j);
                }
            }
}

vector<int> Grid3D::getNearestFluidNeighbor(array<double,3>& x, vector<Particle>& particles, double radius)
{
    Vec position(x[0], x[1], x[2]);
    Vec minBB(x[0]-radius, x[1]-radius, x[2]-radius);
    Vec maxBB(x[0]+radius, x[1]+radius, x[2]+radius);

    array<int,3> minCoord = getGridCoordinate(minBB);
    array<int,3> maxCoord = getGridCoordinate(maxBB);

    minCoord[0] = max(minCoord[0],0);
    minCoord[1] = max(minCoord[1],0);
    minCoord[2] = max(minCoord[2],0);

    maxCoord[0] = min(maxCoord[0], wNumber-1);
    maxCoord[1] = min(maxCoord[1], hNumber-1);
    maxCoord[2] = min(maxCoord[2], dNumber-1);

    vector<int> neighbor;
    for( int _d = minCoord[2]; _d <= maxCoord[2]; ++_d )
    {
        for( int _h = minCoord[1]; _h <= maxCoord[1]; ++_h)
        {
            for( int _w = minCoord[0]; _w <= maxCoord[0]; ++_w)
            {
                int cellId = _w + _h*wNumber+ _d*wNumber*hNumber;
                if( cellId < 0 || cellId >= (int)cells.size() )
                {
                    std::cerr << "Warning : " << cellId << std::endl;
                    neighbor.clear();
                    return neighbor;
                }

                Cell& c = cells[cellId];
                for(int& j : c.fluid)
                {
                    double dist = ( position-particles[j].x ).lengthSquared();
                    if(dist<radius*radius)
                        neighbor.push_back(j);
                }
            }
        }
    }
    return neighbor;
}

bool testGridCreation()
{
    Grid3D grid;
    grid.cellSize = 1;
    Vec min(-1, -1, -1);
    Vec max(1, 1, 1);
    grid.minCorner = min;
    grid.maxCorner = max;
    grid.init();

    bool success = true;

    Vec p;
    int id = -1;

    //Sorting according z position first then y, then x

    p = Vec(-1,-1,-1);
    id = grid.addFluid(0, p);
    if(id != 0){ success = success && false;}

    p = Vec(-1,-1,0);
    id = grid.addFluid(1, p);
    if(id != 1){ success = success && false;}

    p = Vec(-1,0,-1);
    id = grid.addFluid(2, p);
    if(id != 2){ success = success && false;}

    p = Vec(-1,0,0);
    id = grid.addFluid(3, p);
    if(id != 3){ success = success && false;}

    p = Vec(0,-1,-1);
    id = grid.addFluid(4, p);
    if(id != 4){ success = success && false;}

    p = Vec(0,-1,0);
    id = grid.addFluid(5, p);
    if(id != 5){ success = success && false;}

    p = Vec(0,0,-1);
    id = grid.addFluid(6, p);
    if(id != 6){ success = success && false;}

    p = Vec(0,0,0);
    id = grid.addFluid(7, p);
    if(id != 7){ success = success && false;}

    return success;
}

bool testSearchNeighbor()
{
    bool success = true;

    Grid3D grid;
    grid.cellSize = 1;
    Vec min(-2, -2, -2);
    Vec max(2, 2, 2);
    grid.minCorner = min;
    grid.maxCorner = max;
    grid.init();

    Vec target(-2,-2,-2);
    Vec p;
    vector<int> neighbor;

    //Research return the particle itself
    grid.addFluid(0, target);
    neighbor = grid.getNearestFluidNeighbor( target );
    success = success && (neighbor.size() == 1);

    //If particle neighbor is beyond the neighbor cells
    //it is not taken into account
    p = Vec(0,-2,-2);
    grid.addFluid(1,p);

    p = Vec(-2,0,-2);
    grid.addFluid(2,p);

    p = Vec(-2,-2,0);
    grid.addFluid(3,p);

    p = Vec(0,0,0);
    grid.addFluid(4,p);
    neighbor.clear();
    neighbor = grid.getNearestFluidNeighbor( target );
    success = success && (neighbor.size() == 1);

    //As soon as a particle is in the 27 cells arround the cell of the particle
    //It is considered as a neighbor
    p = Vec(-0.1, -2, -2);
    grid.addFluid(5,p);

    p = Vec(-2, -0.1, -2);
    grid.addFluid(6,p);

    p = Vec(-2, -2, -0.1);
    grid.addFluid(7,p);

    p = Vec(-0.1,-0.1,-0.1);
    grid.addFluid(8,p);

    neighbor = grid.getNearestFluidNeighbor( target );
    success = success && (neighbor.size() == 5);

    //Re-run
    grid.clearFluid();
    target = Vec(1,1,1);
    grid.addFluid(0, target);

    //Outside
    p = Vec(-0.1,0,0);
    grid.addFluid(2, p);
    p = Vec(0,-0.1,0);
    grid.addFluid(3, p);
    p = Vec(0,0,-0.1);
    grid.addFluid(4, p);

    neighbor = grid.getNearestFluidNeighbor( target );
    success = success && (neighbor.size() == 1);

    //Inside
    p = Vec(0,0,0);
    grid.addFluid(5, p);
    neighbor = grid.getNearestFluidNeighbor( target );
    success = success && (neighbor.size() == 2);

    //Robustness
    grid.clearFluid();
    grid.addFluid(0, target);
    p = Vec(-2,-2,-2);
    grid.addFluid(1, p);
    p = Vec(-2.1,-2.1,-2.1);
    grid.addFluid(2, p);
    target = Vec(-2.1,-2.1,2.1);
    neighbor = grid.getNearestFluidNeighbor( target );
    success = success && (neighbor.size() == 0);

    return success;
}

void testAll()
{
    int testNumber = 0;
    int testFailure = 0;
    int testSuccess = 0;
    bool success = false;

    success = testGridCreation();
    testNumber++;
    if(success)
    {
        std::cout << "testGridCreation : Success" << std::endl;
        testSuccess++;
    }
    else
    {
        std::cout << "testGridCreation : Failure" << std::endl;
        testFailure++;
    }

    success = testSearchNeighbor();
    testNumber++;
    if(success)
    {
        std::cout << "testSearchNeighbor : Success" << std::endl;
        testSuccess++;
    }
    else
    {
        std::cout << "testSearchNeighbor : Failure" << std::endl;
        testFailure++;
    }

    std::cout << testNumber << " Tests , " << testSuccess << " Success, " << testFailure << " Failure " << std::endl;
}
