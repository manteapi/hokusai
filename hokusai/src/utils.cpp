#include "../include/hokusai/utils.hpp"

AkinciKernel::AkinciKernel()
{
    h = 0;
    m_v1 = 0;
    m_v2 = 0;
}

AkinciKernel::AkinciKernel( double _h)
{
    h = _h;
    m_v1 = 32.0/(M_PI*pow(h,9.0));
    m_v2 = pow(h,6.0)/64.0;
    adhesion = 0.007/pow(h,3.25);
}

AkinciKernel::AkinciKernel( const AkinciKernel& k)
{
    h = k.h;
    m_v1 = k.m_v1;
    m_v2 = k.m_v2;
    adhesion = k.adhesion;
}

AkinciKernel::~AkinciKernel(){}

MonaghanKernel::MonaghanKernel()
{
    h = 0;
    m_v = 0;
    m_g = 0;
}

MonaghanKernel::MonaghanKernel( double _h )
{
    h = _h;
    m_v = 1.0/(4.0*M_PI*h*h*h);
    m_g = 1.0/(4.0*M_PI*h*h*h);
}

MonaghanKernel::MonaghanKernel( const MonaghanKernel& k)
{
    h = k.h;
    m_v = k.m_v;
    m_g = k.m_g;
}

MonaghanKernel::~MonaghanKernel(){}


BoundaryKernel::BoundaryKernel()
{
    h = 0;
    cs = 0;
}

BoundaryKernel::~BoundaryKernel(){}

BoundaryKernel::BoundaryKernel( double _h, double _cs )
{
    h = _h;
    cs = _cs;
}


/*
int mortonNumber( array<int,3>& index )
{
    size_t x = index[0];
    size_t y = index[1];
    size_t z = index[0];
    magnet::math::MortonNumber<3> m(x,y,z);
    size_t morton_num = m.getMortonNum();
    return morton_num;
}
*/

void insertionSort( vector< pair<int,int> >& data )
{
    size_t length = data.size();
    int j = 0;
    pair<int, int> tmp;
    for( size_t i = 0; i < length; ++i )
    {
        j = i;
        while( (j > 0) && (data[j-1].second > data[j].second) )
        {
            tmp = data[j];
            data[j] = data[j-1];
            data[j-1] = tmp;
            j--;
        }
    }
}

Box::Box(){}

Box::Box(Vec& _min, Vec& _max)
{
    min = _min;
    max = _max;
}

Box::~Box(){}

/*void Box::draw()
{
    glLineWidth(3.0);
    glColor4f(0.0,0.0,0.0,1.0);
    glBegin(GL_LINES);
    
    //Square 1
    glVertex3f( min[0], min[1], min[2] ); 
    glVertex3f( max[0], min[1], min[2] ); 
    
    glVertex3f( max[0], min[1], min[2] ); 
    glVertex3f( max[0], max[1], min[2] ); 

    glVertex3f( max[0], max[1], min[2] ); 
    glVertex3f( min[0], max[1], min[2] ); 

    glVertex3f( min[0], max[1], min[2] ); 
    glVertex3f( min[0], min[1], min[2] ); 

    //Square 2
    glVertex3f( min[0], min[1], max[2] ); 
    glVertex3f( max[0], min[1], max[2] ); 
    
    glVertex3f( max[0], min[1], max[2] ); 
    glVertex3f( max[0], max[1], max[2] ); 

    glVertex3f( max[0], max[1], max[2] ); 
    glVertex3f( min[0], max[1], max[2] ); 

    glVertex3f( min[0], max[1], max[2] ); 
    glVertex3f( min[0], min[1], max[2] ); 

    //Square 3
    glVertex3f( min[0], min[1], min[2] ); 
    glVertex3f( min[0], min[1], max[2] ); 
    
    glVertex3f( max[0], min[1], min[2] ); 
    glVertex3f( max[0], min[1], max[2] ); 

    glVertex3f( min[0], max[1], min[2] ); 
    glVertex3f( min[0], max[1], max[2] ); 
    
    glVertex3f( max[0], max[1], min[2] ); 
    glVertex3f( max[0], max[1], max[2] ); 

    glEnd();
}

SolidSphere::SolidSphere(double radius, unsigned int rings, unsigned int sectors)
{
    double const R = 1./(double)(rings-1);
    double const S = 1./(double)(sectors-1);

    vertices.resize(rings * sectors * 3);
    normals.resize(rings * sectors * 3);
    texcoords.resize(rings * sectors * 2);
    std::vector<GLfloat>::iterator v = vertices.begin();
    std::vector<GLfloat>::iterator n = normals.begin();
    std::vector<GLfloat>::iterator t = texcoords.begin();
    for(unsigned int r = 0; r < rings; r++) for(unsigned int s = 0; s < sectors; s++) {
        double const y = sin( -M_PI_2 + M_PI * r * R );
        double const x = cos(2*M_PI * s * S) * sin( M_PI * r * R );
        double const z = sin(2*M_PI * s * S) * sin( M_PI * r * R );

        *t++ = s*S;
        *t++ = r*R;

        *v++ = x * radius;
        *v++ = y * radius;
        *v++ = z * radius;

        *n++ = x;
        *n++ = y;
        *n++ = z;
    }

    indices.resize(rings * sectors * 4);
    std::vector<GLuint>::iterator i = indices.begin();
    for(unsigned int r = 0; r < rings-1; r++) for(unsigned int s = 0; s < sectors-1; s++) {
        *i++ = r * sectors + s;
        *i++ = r * sectors + (s+1);
        *i++ = (r+1) * sectors + (s+1);
        *i++ = (r+1) * sectors + s;
    }
}

SolidSphere::~SolidSphere(){}

void SolidSphere::draw(GLfloat x, GLfloat y, GLfloat z)
{
    glMatrixMode(GL_MODELVIEW);

    glPushMatrix();
    
    glTranslatef(x,y,z);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    //glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    glVertexPointer(3, GL_FLOAT, 0, &vertices[0]);
    glNormalPointer(GL_FLOAT, 0, &normals[0]);
    //glTexCoordPointer(2, GL_FLOAT, 0, &texcoords[0]);
    glDrawElements(GL_QUADS, indices.size(), GL_UNSIGNED_INT, &indices[0]);
    
    glPopMatrix();
}
*/

void write(const char * filename, vector<Vec3<double> > data)
{
    ofstream outputFile;
    outputFile.open(filename);
    outputFile.precision(16);
    for(unsigned int i=0; i <data.size(); ++i)
    {
        outputFile << data[i][0] << " " << data[i][1] <<" " <<  data[i][2] << "\n";
    }
    outputFile.close();
}

static void buildRotationMatrix( float xrad, float yrad, float R[3][3] ) {
    float Rx[3][3] = { {1,0,0},{0,cos(xrad),-sin(xrad)},{0,sin(xrad),cos(xrad)} };
    float Ry[3][3] = { {cos(yrad),0,sin(yrad)}, {0,1,0}, {-sin(yrad),0,cos(yrad)} };
    float Rtmp[3][3];
    for( int i=0; i<3; i++ ) {
        for( int j=0; j<3; j++ ) {
            R[i][j] = Rtmp[i][j] = Rx[i][j];
        }
    }    

    for( int i=0; i<3; i++ ) {
        for( int j=0; j<3; j++ ) {
            R[i][j] = 0.0; 
            for( int k=0; k<3; k++ ) R[i][j] += Rtmp[i][k]*Ry[k][j];
        }
    }    
}

static void transform( float p[3], float R[3][3] ) {
    float p0[3] = { p[0], p[1], p[2] };
    for( int i=0; i<3; i++ ) {
        p[i] = 0.0; 
        for( int k=0; k<3; k++ ) {
            p[i] += R[i][k]*p0[k];
        }
    }    
}



void write_frame(vector<Particle>& particles, int step)
{
    int width = 1024;
    int height = 700;
    float winrate = height/(float)width;

    static float *image = new float[width*height*4];
    static unsigned char *buffer = new unsigned char[width*height*4];
    for( int i=0; i<width; i++ ) {
        for( int j=0; j<height; j++ ) {
            for( int c=0; c<3; c++ ) {
                image[4*(i+j*width)+c] = 0.0; //black
            }
            image[4*(i+j*width)+3] = 1.0; //opacity
        }
    }

    float DENSITY = 0.5;
    int N = 32;
    float blueColor[] = { 0.3, 0.5, 0.8, exp(-DENSITY*N/23.0f) };

    static float R[3][3];
    bool firstTime = true;
    if( firstTime ) {
        //buildRotationMatrix( -0.2, 0.2, R );
        buildRotationMatrix( 0, 0, R );
        firstTime = false;
    }

    // Simple Point-based Rendering
    //float eye = 2.0;
    //float offset = 0.3;
    float eye = 0.02;
    float offset = 0.5;

    for( int n=0; n<particles.size(); n++ ) 
    {
        //if( particles[n]->type == FLUID ) %{
        //float p[3] = { particles[n]->p[0]-0.5, particles[n]->p[1]-0.5, particles[n]->p[2]-0.5 };
        float p[3] = { (float)(particles[n].x[0]-0.5), (float)(particles[n].x[1]-0.5), (float)(particles[n].x[2]-0.5) };
        transform(p,R);
        float z = offset + 0.5 + p[2];
        float x = eye/(eye+z)*(p[0]-0.4);
        float y = eye/(eye+z)*(p[1]+0.25);
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
        //}
    }

    for( int n=0; n<4*width*height; n++ ) {
        buffer[n] = 255.0*fmin(1.0,image[n]);
    }

    char name[256];
    sprintf( name, "frame_%d.bmp", step );
    write_bmp( name, buffer, width, height, true );
}
