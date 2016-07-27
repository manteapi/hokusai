#include "./../include/hokusai/io.hpp"
#include "./../include/hokusai/utils.hpp"
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;
namespace hokusai
{

int countDigits(int number)
{
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

// Skip comments in an ASCII file (lines beginning with #)
void skip_comments(FILE *f)
{
    int c;
    bool in_comment = false;
    while (1) {
        c = fgetc(f);
        if (c == EOF)
            return;
        if (in_comment) {
            if (c == '\n')
                in_comment = false;
        } else if (c == '#') {
            in_comment = true;
        } else if (!isspace(c)) {
            break;
        }
    }
    ungetc(c, f);
}

// Tesselate an arbitrary n-gon.  Appends triangles to "tris".
void tess(const vector<Vec3r > &verts, const vector<int> &thisface,
          vector<Vec3i > &tris)
{
    if (thisface.size() < 3)
        return;
    if (thisface.size() == 3) {
        tris.push_back(Vec3i(thisface[0],
                       thisface[1],
                thisface[2]));
        return;
    }
    if (thisface.size() == 4) {
        // Triangulate in the direction that
        // gives the shorter diagonal
        const Vec3r &p0 = verts[thisface[0]], &p1 = verts[thisface[1]];
        const Vec3r &p2 = verts[thisface[2]], &p3 = verts[thisface[3]];
        HReal d02 = (p0-p2).lengthSquared();
        HReal d13 = (p1-p3).lengthSquared();
        int i = (d02 < d13) ? 0 : 1;
        tris.push_back(Vec3i(thisface[i],
                                  thisface[(i+1)%4],
                       thisface[(i+2)%4]));
        tris.push_back(Vec3i(thisface[i],
                                  thisface[(i+2)%4],
                       thisface[(i+3)%4]));
        return;
    }

    // 5-gon or higher - just tesselate arbitrarily...
    for (size_t i = 2; i < thisface.size(); i++)
        tris.push_back(Vec3i(thisface[0],
                       thisface[i-1],
                thisface[i]));
}

// Read an obj file
bool read_obj(FILE *f, vector< Vec3r >& vertices, vector< Vec3r >& normals, vector< Vec3i >& triangles)
{


    vector<int> thisface;
    while (1) {
        skip_comments(f);
        if (feof(f))
            return true;
        char buf[1024];
        GET_LINE();
        if (LINE_IS("v ") || LINE_IS("v\t")) {
            HReal x, y, z;
            if (sscanf(buf+1, "%lf %lf %lf", &x, &y, &z) != 3) {
                return false;
            }
            vertices.push_back(Vec3r(x,y,z));
        } else if (LINE_IS("vn ") || LINE_IS("vn\t")) {
            HReal x, y, z;
            if (sscanf(buf+2, "%lf %lf %lf", &x, &y, &z) != 3) {
                return false;
            }
            normals.push_back(Vec3r(x,y,z));
        } else if (LINE_IS("f ") || LINE_IS("f\t") ||
                   LINE_IS("t ") || LINE_IS("t\t")) {
            thisface.clear();
            char *c = buf;
            while (1) {
                while (*c && *c != '\n' && !isspace(*c))
                    c++;
                while (*c && isspace(*c))
                    c++;
                int thisf;
                if (sscanf(c, " %d", &thisf) != 1)
                    break;
                if (thisf < 0)
                    thisf += vertices.size();
                else
                    thisf--;
                thisface.push_back(thisf);
            }
            tess(vertices, thisface, triangles);
        }
    }

    // XXX - FIXME
    // Right now, handling of normals is fragile: we assume that
    // if we have the same number of normals as vertices,
    // the file just uses per-vertex normals.  Otherwise, we can't
    // handle it.
    if (vertices.size() != normals.size())
        normals.clear();

    return true;
}

HokusaiImporter::~HokusaiImporter(){}

HokusaiImporter::HokusaiImporter(const std::string& fileName, System& system)
{
    std::ifstream file(fileName.c_str(), ios::out|ios::binary);
    if(file)
    {
        //Particle number
        int particleNumber;
        file.read((char*)&particleNumber, sizeof(particleNumber));
        //Smoothing radius
        HReal smoothingRadius;
        file.read((char*)&smoothingRadius, sizeof(smoothingRadius));
        //Particle positions
        float x,y,z;
        for(int i=0; i<particleNumber; ++i)
        {
            file.read((char*)&x, sizeof(x));
            file.read((char*)&y, sizeof(y));
            file.read((char*)&z, sizeof(z));
            system.addFluidParticle(Vec3r(x,y,z), Vec3r(0,0,0));
        }
        file.close();
    }
    else
    {
        std::cout << "HokusaiExporter:: error while reading the file" << std::endl;
    }
}

HokusaiExporter::~HokusaiExporter(){}

HokusaiExporter::HokusaiExporter(const std::string& fileName, const System& system)
{
    std::ofstream file(fileName.c_str(), ios::out|ios::binary);
    if(file)
    {
        //Particle number
        const int particleNumber = system.particleNumber();
        file.write((char*)&particleNumber, sizeof(particleNumber));
        //Smoothing radius
        const HReal smoothingRadius = system.fluidParams().smoothingRadius();
        file.write((char*)&smoothingRadius, sizeof(smoothingRadius));
        //Particle positions
        const std::vector<Particle>& particles = system.particles();
        for(int i=0; i<particleNumber; ++i)
        {
            const Particle& p = particles[i];
            float x=p.x[0], y=p.x[1], z=p.x[2];
            file.write((char*)&x, sizeof(x));
            file.write((char*)&y, sizeof(y));
            file.write((char*)&z, sizeof(z));
        }
        file.close();
    }
    else
    {
        std::cout << "HokusaiExporter:: error while writting the file" << std::endl;
    }
}

BlenderExporter::~BlenderExporter()
{
}

BlenderExporter::BlenderExporter()
{
}

BlenderExporter::BlenderExporter(const std::string& prefix, const std::string& path, 
                                 const int& particleNumber, const int& frameNumber)
{
    m_particleNumber = particleNumber;
    m_prefix = prefix;
    m_path = path;
    m_frameNumber = frameNumber;
    m_frameCounter = 0;

    std::stringstream padding; 
    padding.fill('0'); 
    padding.width( max(6,countDigits(particleNumber)) );
    padding << '0';
    std::string filename = m_path + m_prefix + "_" + padding.str() + "_00.bphys";

    std::ofstream file(filename.c_str(), ios::out|ios::binary);
    if(file)
    {
        //HEADER
        const std::string bphysics = "BPHYSICS";
        file.write(bphysics.c_str(), sizeof(bphysics));
        //FLAG (compression algorithm maybe)
        unsigned int flag = 1;
        file.write((char*)&flag, sizeof(flag));
        //NUMBER OF DATA
        unsigned int size = particleNumber;
        file.write((char*)&size, sizeof(size));
        //dunno
        unsigned int unknown = 64;
        file.write((char*)&unknown, sizeof(unknown));
        //DATA
        for(int i=0; i<particleNumber; ++i)
        {
            float frameStart = 1;
            float frameEnd = m_frameNumber+1;
            float frameNumber = m_frameNumber;
            file.write((char*)&frameStart, sizeof(frameStart));
            file.write((char*)&frameEnd, sizeof(frameEnd));
            file.write((char*)&frameNumber, sizeof(frameNumber));
        }
        file.close();
    }
}

void BlenderExporter::debugHeader(const std::string& fileName)
{
    std::ifstream readStream(fileName.c_str(), ios::binary | ios::in);
    if(readStream)
    {
        char bphysics[8];
        readStream.read((char*)&bphysics, sizeof(bphysics));
        std::cout << "Header: " << bphysics << std::endl;

        int flag;
        readStream.read((char*)&flag, sizeof(flag));
        std::cout << "Flag: " << flag << std::endl;

        int particleNumber;
        readStream.read((char*)&particleNumber, sizeof(particleNumber));
        std::cout << "ParticleNumber: " << particleNumber << std::endl;

        int unknown_64;
        readStream.read((char*)&unknown_64, sizeof(unknown_64));
        std::cout << "Unknown_64: " << unknown_64 << std::endl;

        for(int i=0; i<particleNumber; ++i)
        {
            float frameStart=0.0, frameEnd=0.0, frameNumber=0.0;
            readStream.read((char*)&frameStart, sizeof(frameStart));
            readStream.read((char*)&frameEnd, sizeof(frameEnd));
            readStream.read((char*)&frameNumber, sizeof(frameNumber));
            std::cout << "Particle " << i << ", lifeTime :" << frameStart << ", " << frameEnd << ", " << frameNumber << std::endl;
        }

        readStream.close();
    }
    else
    {
        std::cerr << "Reading : Unable to open the file" << std::endl;
        std::cerr << "Filename :" << fileName << std::endl;
    }
}

void BlenderExporter::debugFrame(const std::string& fileName)
{
    std::ifstream readStream(fileName.c_str(), ios::binary | ios::in);
    if(readStream)
    {
        char bphysics[8];
        readStream.read((char*)&bphysics, sizeof(bphysics));
        std::cout << "Header: " << bphysics << std::endl;

        unsigned int type;
        readStream.read((char*)&type, sizeof(type));
        std::cout << "Type: " << type << std::endl;

        unsigned int particleNumber;
        readStream.read((char*)&particleNumber, sizeof(particleNumber));
        std::cout << "ParticleNumber: " << particleNumber << std::endl;

        unsigned int dataType;
        readStream.read((char*)&dataType, sizeof(dataType));
        std::cout << "DataType: " << dataType << std::endl;

        for(unsigned int i=0; i<particleNumber; ++i)
        {
            float x1, x2, x3, v1, v2, v3;
            readStream.read((char*)&x1, sizeof(x1));
            readStream.read((char*)&x2, sizeof(x2));
            readStream.read((char*)&x3, sizeof(x3));
            readStream.read((char*)&v1, sizeof(v1));
            readStream.read((char*)&v2, sizeof(v2));
            readStream.read((char*)&v3, sizeof(v3));
            std::cout << "Particle " << i << " : x(" << x1 << ", " << x2 << ", " << x3 << "), v(" << v1 << ", " << v2 << ", " << v3 << ")" << std::endl;
        }

        readStream.close();
    }
    else
    {
        std::cerr << "Reading : Unable to open the file" << std::endl;
        std::cerr << "Filename :" << fileName << std::endl;
    }
}

void BlenderExporter::apply(const System& system)
{
    //Build a file per frame
    std::stringstream padding; 
    padding.fill('0'); 
    padding.width( max(6,countDigits(m_particleNumber)) );
    m_frameCounter++;
    padding << std::to_string(m_frameCounter);
    std::string filename = m_path + m_prefix + "_" + padding.str() + "_00.bphys";
    std::ofstream file(filename.c_str(), ios::out|ios::binary);
    if(file)
    {
        //HEADER
        const string bphysics = "BPHYSICS";
        file.write(bphysics.c_str(), sizeof(bphysics));
        //TYPE
        unsigned int type = 1;
        file.write((char*)&type,sizeof(type));
        //NUMBER OF DATA
        unsigned int size = m_particleNumber;
        file.write((char*)&size,sizeof(size));
        //DATA TYPE (index + location + velocity = 7)
        unsigned int dataType = 7;
        file.write((char*)&dataType,sizeof(dataType));
        for(int i=0; i<system.particleNumber(); ++i)
        {
            const Particle& pi = system.particles()[i];
            float x1 = pi.x[0], x2 = pi.x[1], x3 = pi.x[2];
            float v1 = pi.v[0], v2 = pi.v[1], v3 = pi.v[2];
            //DATA
            //Index
            file.write((char*)&i,sizeof(i));
            //Location
            file.write((char*)&x1, sizeof(x1));
            file.write((char*)&x2, sizeof(x2));
            file.write((char*)&x3, sizeof(x3));
            //Velocity
            file.write((char*)&v1, sizeof(v1));
            file.write((char*)&v2, sizeof(v2));
            file.write((char*)&v3, sizeof(v3));
            //Orientation
            //file.write((char*)orientationblender, 16);
        }
        file.close();
    }
}

std::string& BlenderExporter::prefix()
{
    return m_prefix;
}

const std::string& BlenderExporter::prefix() const
{
    return m_prefix;
}

std::string& BlenderExporter::path()
{
    return m_path;
}

const std::string& BlenderExporter::path() const
{
    return m_path;
}

int& BlenderExporter::frameNumber()
{
    return m_frameNumber;
}

const int& BlenderExporter::frameNumber() const
{
    return m_frameNumber;
}

int& BlenderExporter::particleNumber()
{
    return m_particleNumber;
}

const int& BlenderExporter::particleNumber() const
{
    return m_particleNumber;
}

void write(const char * filename, std::vector<HReal> data)
{
    std::ofstream outputFile;
    outputFile.open(filename);
    outputFile.precision(16);
    for(unsigned int i=0; i <data.size(); ++i)
    {
        outputFile << data[i] <<"\n";
    }
    outputFile.close();
}

void write(const char * filename, std::vector<Vec3r > data)
{
    std::ofstream outputFile;
    outputFile.open(filename);
    outputFile.precision(16);
    for(unsigned int i=0; i <data.size(); ++i)
    {
        outputFile << data[i][0] << " " << data[i][1] <<" " <<  data[i][2] << "\n";
    }
    outputFile.close();
}

void write_frame(const System& sph, int step, HReal offset)
{
    int width = 1024;
    int height = 700;
    HReal winrate = height/(HReal)width;

    static HReal *image = new HReal[width*height*4];
    static unsigned char *buffer = new unsigned char[width*height*4];
    for( int i=0; i<width; i++ ) {
        for( int j=0; j<height; j++ ) {
            for( int c=0; c<3; c++ ) {
                image[4*(i+j*width)+c] = 0.0; //black
            }
            image[4*(i+j*width)+3] = 1.0; //opacity
        }
    }

    HReal DENSITY = 0.5;
    int N = 32;
    HReal blueColor[] = { 0.3, 0.5, 0.8, exp(-DENSITY*N/23.0f) };
    HReal yellowColor[] = { 0.8, 0.2, 0.0, exp(-DENSITY*N/23.0f) };

    static HReal R[3][3];
    bool firstTime = true;
    if( firstTime ) {
        //buildRotationMatrix( -0.2, 0.2, R );
        buildRotationMatrix( 0, 0, R );
        firstTime = false;
    }

    // Simple Point-based Rendering
    //HReal eye = 2.0;
    //HReal offset = 0.3;
    HReal eye = 2.0;

    const std::vector<Particle>& particles = sph.particles();
    for( size_t n=0; n<particles.size(); n++ )
    {
        HReal p[3] = { (HReal)(particles[n].x[0]-0.5), (HReal)(particles[n].x[1]-0.5), (HReal)(particles[n].x[2]-0.5) };
        transform(p,R);
        HReal z = offset + 0.5 + p[2];
        HReal x = eye/(eye+z)*(p[0]-0.4);
        HReal y = eye/(eye+z)*(p[1]+0.25);
        if( x > -0.5 && x < 0.5 && y > -winrate/2.0 && y < winrate/2.0 ) {
            int i = width*(x+0.5);
            int j = width*(y+winrate/2.0);
            int w = 1;
            for( int ni=i-w; ni<i+1; ni++ ) {
                for( int nj=j-w; nj<j+1; nj++ ) {
                    if( ni>0 && ni<width && nj>0 && nj<height ) {
                        for( int c=0; c<3; c++ ) {
                            image[4*(ni+nj*width)+c] = blueColor[3]*blueColor[c] + (1.0-blueColor[3])*image[4*(ni+nj*width)+c];
                        }
                    }
                }
            }
        }
    }

    const std::vector<Boundary>& boundaries= sph.boundaries();
    for( size_t n=0; n<boundaries.size(); n++ )
    {
        HReal p[3] = { (HReal)(boundaries[n].x[0]-0.5), (HReal)(boundaries[n].x[1]-0.5), (HReal)(boundaries[n].x[2]-0.5) };
        transform(p,R);
        HReal z = offset + 0.5 + p[2];
        HReal x = eye/(eye+z)*(p[0]-0.4);
        HReal y = eye/(eye+z)*(p[1]+0.25);
        if( x > -0.5 && x < 0.5 && y > -winrate/2.0 && y < winrate/2.0 ) {
            int i = width*(x+0.5);
            int j = width*(y+winrate/2.0);
            int w = 1;
            for( int ni=i-w; ni<i+1; ni++ ) {
                for( int nj=j-w; nj<j+1; nj++ ) {
                    if( ni>0 && ni<width && nj>0 && nj<height ) {
                        for( int c=0; c<3; c++ ) {
                            image[4*(ni+nj*width)+c] = yellowColor[3]*yellowColor[c] + (1.0-yellowColor[3])*image[4*(ni+nj*width)+c];
                        }
                    }
                }
            }
        }
    }

    for( int n=0; n<4*width*height; n++ ) {
        buffer[n] = 255.0*fmin(1.0,image[n]);
    }

    char name[256];
    sprintf( name, "frame_%d.bmp", step );
    write_bmp( name, buffer, width, height, true );
}

}//namespace hokusai
