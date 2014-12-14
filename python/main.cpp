#include <limits>
#include <cmath>
#include <array>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string.h>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

int countDigits(int number) {
    if (number < 10) {
        return 1;
    }
    int count = 0;
    while (number > 0) {
        number /= 10;
        count++;
    }
    return count;
}

void toBlenderVector(float* dest, Vector3d& src)
{
    dest[0] = (float)src[0];
    dest[1] = (float)src[1];
    dest[2] = (float)src[2];
}

void writeBphysHeader(const string& base, int particleNumber, int beginFrame, int endFrame, int lifeTime)
{
    //Build a file per frame
    stringstream padding; padding.fill('0'); padding.width( max(6,countDigits(particleNumber)) );
    padding << '0';
    string filename = base + "_" + padding.str() + "_00.bphys";
    ofstream file(filename.c_str(), ios::out|ios::binary);
    if(file)
    {
        //HEADER
        const string bphysics = "BPHYSICS";
        file.write(bphysics.c_str(), 8);
        //FLAG (compression algorithm maybe)
        unsigned int flag = 1;
        file.write((char*)&flag,4);
        //NUMBER OF DATA
        unsigned int size = particleNumber;
        file.write((char*)&size,4);
        //dunno
        unsigned int unknown = 64;
        file.write((char*)&unknown,4);
        //DATA
        for(int i=0; i<particleNumber; ++i)
        {
            float frame_info[3] = {(float)beginFrame, (float)endFrame, (float)lifeTime};
            file.write((char*)&frame_info, 12);
        }
        file.close();
    }
}

void writeBphysData(string& base, int frameNumber, vector<Vector3d>& position, vector<Vector3d>& velocity)
{
    //Build a file per frame
    stringstream padding; padding.fill('0'); padding.width( max(6,countDigits(position.size())) );
    padding << (frameNumber+1);
    string filename = base + "_" + padding.str() + "_00.bphys";
    ofstream file(filename.c_str(), ios::out|ios::binary);
    if(file)
    {
        //HEADER
        const string bphysics = "BPHYSICS";
        file.write(bphysics.c_str(), 8);
        //TYPE
        unsigned int type = 1;
        file.write((char*)&type,4);
        //NUMBER OF DATA
        unsigned int size = position.size();
        file.write((char*)&size,4);
        //DATA TYPE (index + location + velocity = 7)
        unsigned int dataType = 7;
        file.write((char*)&dataType,4);
        for(unsigned int i=0; i<position.size(); ++i)
        {
            //DATA
            //float frame = (float) frameNumber;
            float posblender[3], velblender[3];
            //float orientationblender[4] = {0,0,0,0};
            toBlenderVector(posblender, position[i]);
            toBlenderVector(velblender, velocity[i]);
            //Index
            file.write((char*)&i,4);
            //Location
            file.write((char*)posblender, 12);
            //Velocity
            file.write((char*)velblender, 12);
            //Orientation
            //file.write((char*)orientationblender, 16);
        }
        file.close();
    }
}

//Required information in the .blend :
//  particle number
//  starting frame, ending frame
int main()
{
    string base = "toto";
    int frameNumber = 0;
    int beginFrame = 1, endFrame = 11, lifeTime = 10;
    vector<Vector3d> position, velocity;
    position.push_back(Vector3d(-1,1,0));
    position.push_back(Vector3d(1,1,0));
    velocity.push_back(Vector3d(0,1,0));
    velocity.push_back(Vector3d(0,1,0));

    writeBphysHeader(base, position.size(), beginFrame, endFrame, lifeTime);
    for(int i=0; i<lifeTime; ++i)
    {
        writeBphysData(base,frameNumber, position, velocity);
        frameNumber++;
        for(unsigned int j=0; j<position.size(); ++j)
        {
            position[j] = position[j] + ( 1.0/24.0)*velocity[j];
        }
    }
}
