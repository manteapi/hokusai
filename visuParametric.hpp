#ifndef VISUPARAMETRIC_HPP
#define VISUPARAMETRIC_HPP

#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/type_precision.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/norm.hpp>
#include <GL/glew.h>
#include <array>
#include "Vec.hpp"
#include "renderable.hpp"

using namespace std;
using namespace glm;

typedef glm::detail::tvec3<int> vec3i;

class VisuParametric : public Renderable
{
    public :
    VisuParametric();
    ~VisuParametric();

    virtual void initShader( int shaderId );
    virtual void createvbo();
    virtual void reloadvbo();
    void drawVisuParametric();
    void drawVoxel();
    void rasterizeLinesExample();
    void rasterizeLineExample();
    void rasterizeTriangleExample1();
    void rasterizeTriangleExample2();
    virtual void drawShader();
    virtual void drawFix();
    virtual float rayPicking(vec3& source, vec3& ray);
    virtual const mat4& getModel() { return model; }
    virtual const GLuint getModelLoc() { return modelLoc; }
    virtual void keyPressedEvent(sf::Event& e);
    void getCubeData(std::vector<vec3>& points, std::vector<vec4>& colors, const vec3& defaultPosition=vec3(0,0,0), const float& scale=1.0, const vec4& defaultColor=vec4(-1,-1,-1,-1));

    std::vector<Vec3<double> > getDiskSampling(const Vec3<double>& center, double radius, double spacingX, double spacingY);
    std::vector<Vec3<double> > getEllipsoidSampling(const Vec3<double>& center, double axis_1, double axis_2, double axis_3, double spacingX, double spacingY);
    std::vector<Vec3<double> > getSphereSampling(const Vec3<double>& center, double radius, double spacingX, double spacingY);
    std::vector<Vec3<double> > getHemiSphereSampling(const Vec3<double>& center, double radius, double spacingX, double spacingY);
    std::vector<Vec3<double> > getCapsuleSampling(const Vec3<double>& center, double radius, double height, double spacingX, double spacingY);
    std::vector< Vec3<double> > getTorusSampling(const Vec3<double>& center, double tubeRadius, double innerRadius, double spacingX, double spacingY);
    std::vector< Vec3<double> > getConeSampling(const Vec3<double>& center, double height, double stopHeight, double baseRadius, double spacingX, double spacingY);
    std::vector< Vec3<double> > getCylinderSampling(const Vec3<double>& center, double height, double baseRadius, double spacingX, double spacingY);
    std::vector<Vec3<double> > getDiskSampling(const Vec3<double>& center, double radius, double spacing);
    std::vector<Vec3<double> > getPyramidSampling(const Vec3<double>& _center, double base, double height, double spacing);
private:
    mat4 model;
    vec3 center;
    float radius;
    int picked;

    public :
    //Data
    
    vector<vec3> vertices;
    vector<vec4> colors;

    //Buffer id

    GLuint positionBuffer;
    GLuint colorBuffer;

    //Location

    GLuint modelLoc;
    GLuint positionLoc;
    GLuint colorLoc;

    vec3 getVec(const Vec3<double>& v)
    {
        vec3 r;
        for(int i=0; i<3; ++i)
            r[i] = v[i];
        return r;
    }

    vec3 getVec(const Vec3<int>& v)
    {
        vec3 r;
        for(int i=0; i<3; ++i)
            r[i] = v[i];
        return r;
    }
};
#endif
