#include <iostream>
#include <string>

#include <GL/glew.h>
#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>
#include <SFML/Graphics/Font.hpp>
#include <SFML/System/String.hpp>
#include <glm/gtc/type_precision.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "shader.hpp"

#include <hokusai/utility.hpp>
#include <hokusai/system.hpp>
#include <hokusai/implicitFunction.hpp>
#include <hokusai/gridUtility.hpp>
#include <hokusai/marchingSquare.hpp>
#include <hokusai/io.hpp>
#include <AntTweakBar.h>

using namespace std;

void TW_CALL cleanParticleCB(void * clientData)
{
    hokusai::System* tmp = static_cast<hokusai::System*>(clientData);
    tmp->cleanFluidParticle();
}

void TW_CALL damBreakCB(void * clientData)
{
    hokusai::System* tmp = static_cast<hokusai::System*>(clientData);
    Vec2d fluidOffset(-0.25,-0.25);
    Vec2d fluidBox(0.5,0.5);
    tmp->addParticleBox(fluidOffset, fluidBox);
}

void TW_CALL updateParamCB(void * clientData)
{
    hokusai::System* tmp = static_cast<hokusai::System*>(clientData);
    tmp->updateParameters();
}

int main()
{
    //Set simulation
    int resolution = 2e3;
    hokusai::System sph(resolution);

    ///Set fluid
    Vec2d fluidOffset(-0.25,-0.45);
    Vec2d fluidBox(0.5,0.5);
    sph.addParticleBox(fluidOffset, fluidBox);

    //Set boundary
    Vec2d securityOffset(1.05*sph.getSmoothingRadius());
    Vec2d boundOffset(-0.5,-0.5);
    boundOffset -= 1.0*securityOffset;
    Vec2d boundBox(1.0,1.0);
    boundBox += securityOffset;
    sph.addBoundaryBox(boundOffset, boundBox);

    ///Initialize parameters
    sph.init();
    sph.setTimeStep(2e-3);
    sph.setFluidCohesion(0);
    sph.setViscosity(1e-3);
    sph.setBoundaryAdhesion(0);
    sph.setBoundaryFriction(0);
//    sph.setGravity(Vec2d(0.0,0.0));
    sph.debugFluid();
    sph.gridInfo.info();

    std::vector<double> scalarField;
    hokusai::Grid2dUtility gridInfo;
    double surfaceResolution = 5e-3;
    double initialValue = 1.0;
    double isovalue = 0.0;
    std::vector<hokusai::Edge> edges;
    hokusai::computeScalarField(scalarField, gridInfo, sph, surfaceResolution, initialValue);
    hokusai::polygonize(edges, scalarField, gridInfo, isovalue);

    //Initialize visualization data
    std::vector< glm::vec2 > particle_center;
    std::vector< glm::vec4 > particle_color;
    std::vector< Vec2d > tmp_particle_center;
    tmp_particle_center = sph.getFluidPosition();
    for(size_t i=0; i<tmp_particle_center.size(); ++i)
    {
        particle_center.push_back(glm::vec2(tmp_particle_center[i][0], tmp_particle_center[i][1]));
        particle_color.push_back(glm::vec4(0.0,0.0,1.0,1.0));
    }
    tmp_particle_center = sph.getBoundaryPosition();
    for(size_t i=0; i<tmp_particle_center.size(); ++i)
    {
        particle_center.push_back(glm::vec2(tmp_particle_center[i][0], tmp_particle_center[i][1]));
        particle_color.push_back(glm::vec4(1.0,0.0,0.0,1.0));
    }

    std::vector< glm::vec2 > grid_vertices;
    std::vector< glm::vec4 > grid_color;
    Vec2d gOffset = sph.gridInfo.offset;
    Vec2d gScale = sph.gridInfo.scale;
    grid_vertices.push_back( glm::vec2(gOffset[0], gOffset[1]) );
    grid_vertices.push_back( glm::vec2(gOffset[0]+gScale[0], gOffset[1]) );
    grid_vertices.push_back( glm::vec2(gOffset[0]+gScale[0], gOffset[1]+gScale[1]) );
    grid_vertices.push_back( glm::vec2(gOffset[0], gOffset[1]+gScale[1]) );
    grid_color.push_back( glm::vec4(0.0,0.0,0.0,1.0) );
    grid_color.push_back( glm::vec4(0.0,0.0,0.0,1.0) );
    grid_color.push_back( glm::vec4(0.0,0.0,0.0,1.0) );
    grid_color.push_back( glm::vec4(0.0,0.0,0.0,1.0) );

    std::vector< glm::vec2 > fluidSurface_vertices;
    std::vector< glm::vec4 > fluidSurface_color;
    for(size_t i=0; i<edges.size(); ++i)
    {
        fluidSurface_vertices.push_back( glm::vec2(edges[i].p1[0], edges[i].p1[1]) );
        fluidSurface_vertices.push_back( glm::vec2(edges[i].p2[0], edges[i].p2[1]) );
        fluidSurface_color.push_back( glm::vec4(1.0,1.0,0.0,1.0) );
        fluidSurface_color.push_back( glm::vec4(1.0,1.0,0.0,1.0) );
    }

    ///Window creation
    sf::RenderWindow window(sf::VideoMode(1280,720), "fluid", sf::Style::Default, sf::ContextSettings(32));

    ///OpenGL context
    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if( GLEW_OK != err )
    {
        cout << "Error : " << glewGetErrorString(err) << endl;
    }
    cout << "Status : Using GLEW " << glewGetString(GLEW_VERSION) << std::endl;

    ///OpenGL info
    const GLubyte* renderer = glGetString (GL_RENDERER); // get renderer string
    const GLubyte* version = glGetString (GL_VERSION); // version as a string
    std::cout << "Renderer : " << renderer << std::endl;
    std::cout << "OpenGL version supported : " << version << std::endl;

    TwInit(TW_OPENGL, NULL);
    TwWindowSize(window.getSize().x, window.getSize().y);
    TwBar *Tbar = TwNewBar("Tweak");
    TwAddVarRW(Tbar, "surf. res.", TW_TYPE_DOUBLE, &surfaceResolution, "min=1e-4 max=1e-0 step=1e-4 group=Surface label='surface resolution' ");
    TwAddVarRW(Tbar, "viscosity", TW_TYPE_DOUBLE, &sph.getViscosity(), "min=0.0 max=1e1 step=1e-3 group=Fluid label='viscosity' ");
    TwAddVarRW(Tbar, "cohesion", TW_TYPE_DOUBLE, &sph.getFluidCohesion(), "min=0.0 max=1e1 step=1e-3 group=Fluid label='cohesion' ");
    TwAddVarRW(Tbar, "restDensity", TW_TYPE_DOUBLE, &sph.getRestDensity(), "min=0.0 max=1e5 step=1e1 group=Fluid label='restDensity' ");
    TwAddButton(Tbar, "cleanParticle", cleanParticleCB, &sph, " group=Fluid label='clean fluid' ");
    TwAddButton(Tbar, "damBreak", damBreakCB, &sph, " group=Fluid label='dam break' ");
    TwAddButton(Tbar, "updateParam", updateParamCB, &sph, " group=Fluid label='update param.' ");
    TwAddVarRW(Tbar, "adhesion", TW_TYPE_DOUBLE, &sph.getBoundaryAdhesion(), "min=0.0 max=1e1 step=1e-3 group=Boundary label='adhesion' ");
    TwAddVarRW(Tbar, "friction", TW_TYPE_DOUBLE, &sph.getBoundaryFriction(), "min=0.0 max=1e1 step=1e-3 group=Boundary label='friction' ");
    TwAddVarRO(Tbar, "volume", TW_TYPE_DOUBLE, &sph.getRealVolumeValue(), "group=Stats label='volume' ");
    TwAddVarRO(Tbar, "density", TW_TYPE_DOUBLE, &sph.getMeanDensityValue(), "group=Stats label='density' ");
    TwAddVarRO(Tbar, "densityFluct.", TW_TYPE_DOUBLE, &sph.getDensityFluctuationValue(), "group=Stats label='densityFluct.' ");

    //glm::mat4 projMat = glm::mat4(1.0f);
    //glm::mat4 projMat = glm::ortho(0.0f,(float)window.getSize().x,0.0f,(float)window.getSize().y);
    double windowRatio = window.getSize().x/(double)window.getSize().y;
    glm::mat4 projMat = glm::ortho((float)(-1.0f*windowRatio),(float)(1.0f*windowRatio),-1.0f,1.0f);

    //Global settings
    glClearColor(0.0f,0.5f,0.5f,1.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_POINT_SPRITE);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    //Compile Shader and link shader
    Shader particleShader;
    string particleVertexPath = "./../particle_vertex.glsl";
    string particleFragmentPath = "./../particle_fragment.glsl";
    particleShader.load(particleVertexPath.c_str(), particleFragmentPath.c_str());

    Shader gridShader;
    string gridVertexPath = "./../grid_vertex.glsl";
    string gridFragmentPath = "./../grid_fragment.glsl";
    gridShader.load(gridVertexPath.c_str(), gridFragmentPath.c_str());

    Shader fluidSurfaceShader;
    string fluidSurfaceVertexPath = "./../fluidSurface_vertex.glsl";
    string fluidSurfaceFragmentPath = "./../fluidSurface_fragment.glsl";
    fluidSurfaceShader.load(fluidSurfaceVertexPath.c_str(), fluidSurfaceFragmentPath.c_str());

    //Set Shader variables
    GLuint particleProjMatLoc = glGetUniformLocation(particleShader.id(), "projMat");
    GLuint gridProjMatLoc = glGetUniformLocation(gridShader.id(), "projMat");
    GLuint fluidSurfaceProjMatLoc = glGetUniformLocation(fluidSurfaceShader.id(), "projMat");

    GLuint particlePositionLoc = glGetAttribLocation(particleShader.id(), "particle_position");
    GLuint particleColorLoc = glGetAttribLocation(particleShader.id(), "particle_color");
    GLuint particlePositionBuffer;
    glGenBuffers(1, &particlePositionBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, particlePositionBuffer);
    glBufferData(GL_ARRAY_BUFFER, particle_center.size()*sizeof(glm::vec2), particle_center.data(), GL_STATIC_DRAW);
    GLuint particleColorBuffer;
    glGenBuffers(1, &particleColorBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, particleColorBuffer);
    glBufferData(GL_ARRAY_BUFFER, particle_color.size()*sizeof(glm::vec4), particle_color.data(), GL_STATIC_DRAW);

    GLuint gridPositionLoc = glGetAttribLocation(gridShader.id(), "grid_position");
    GLuint gridColorLoc = glGetAttribLocation(gridShader.id(), "grid_color");
    GLuint gridPositionBuffer;
    glGenBuffers(1, &gridPositionBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, gridPositionBuffer);
    glBufferData(GL_ARRAY_BUFFER, grid_vertices.size()*sizeof(glm::vec2), grid_vertices.data(), GL_STATIC_DRAW);
    GLuint gridColorBuffer;
    glGenBuffers(1, &gridColorBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, gridColorBuffer);
    glBufferData(GL_ARRAY_BUFFER, grid_color.size()*sizeof(glm::vec4), grid_color.data(), GL_STATIC_DRAW);

    GLuint fluidSurfacePositionLoc = glGetAttribLocation(fluidSurfaceShader.id(), "fluidSurface_position");
    GLuint fluidSurfaceColorLoc = glGetAttribLocation(fluidSurfaceShader.id(), "fluidSurface_color");
    GLuint fluidSurfacePositionBuffer;
    glGenBuffers(1, &fluidSurfacePositionBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, fluidSurfacePositionBuffer);
    glBufferData(GL_ARRAY_BUFFER, fluidSurface_vertices.size()*sizeof(glm::vec2), fluidSurface_vertices.data(), GL_STATIC_DRAW);
    GLuint fluidSurfaceColorBuffer;
    glGenBuffers(1, &fluidSurfaceColorBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, fluidSurfaceColorBuffer);
    glBufferData(GL_ARRAY_BUFFER, fluidSurface_color.size()*sizeof(glm::vec4), fluidSurface_color.data(), GL_STATIC_DRAW);

    bool run = true;
    while(run)
    {
        //Simulate
        sph.simulate();

        //Gather visualization data
        particle_center.clear();
        particle_color.clear();
        tmp_particle_center = sph.getFluidPosition();
        for(size_t i=0; i<tmp_particle_center.size(); ++i)
        {
            particle_center.push_back(glm::vec2(tmp_particle_center[i][0], tmp_particle_center[i][1]));
            particle_color.push_back(glm::vec4(0.0,0.0,1.0,1.0));
        }
        tmp_particle_center = sph.getBoundaryPosition();
        for(size_t i=0; i<tmp_particle_center.size(); ++i)
        {
            particle_center.push_back(glm::vec2(tmp_particle_center[i][0], tmp_particle_center[i][1]));
            particle_color.push_back(glm::vec4(1.0,0.0,0.0,1.0));
        }

        edges.clear();
        scalarField.clear();
        hokusai::computeScalarField(scalarField, gridInfo, sph, surfaceResolution, initialValue);
        hokusai::polygonize(edges, scalarField, gridInfo, isovalue);
        //hokusai::exportOBJ("toto.obj", edges);

        fluidSurface_vertices.clear();
        fluidSurface_color.clear();
        for(size_t i=0; i<edges.size(); ++i)
        {
            fluidSurface_vertices.push_back( glm::vec2(edges[i].p1[0], edges[i].p1[1]) );
            fluidSurface_vertices.push_back( glm::vec2(edges[i].p2[0], edges[i].p2[1]) );
            fluidSurface_color.push_back( glm::vec4(1.0,1.0,0.0,1.0) );
            fluidSurface_color.push_back( glm::vec4(1.0,1.0,0.0,1.0) );
        }

        //Reload vbo
        glUseProgram(particleShader.id());
        glBindBuffer(GL_ARRAY_BUFFER, particlePositionBuffer);
        glBufferData(GL_ARRAY_BUFFER, particle_center.size()*sizeof(glm::vec2), particle_center.data(), GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, particleColorBuffer);
        glBufferData(GL_ARRAY_BUFFER, particle_color.size()*sizeof(glm::vec4), particle_color.data(), GL_STATIC_DRAW);
        glUseProgram(0);

        glUseProgram(fluidSurfaceShader.id());
        glBindBuffer(GL_ARRAY_BUFFER, fluidSurfacePositionBuffer);
        glBufferData(GL_ARRAY_BUFFER, fluidSurface_vertices.size()*sizeof(glm::vec2), fluidSurface_vertices.data(), GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, fluidSurfaceColorBuffer);
        glBufferData(GL_ARRAY_BUFFER, fluidSurface_color.size()*sizeof(glm::vec4), fluidSurface_color.data(), GL_STATIC_DRAW);
        glUseProgram(0);

        ///Clear screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        ///Catch event
        sf::Event event;
        while(window.pollEvent(event))
        {
            // send event to AntTweakBar
            int handled = TwEventSFML(&event, 1, 6); // assume SFML version 1.6 here
            if( !handled ) // if event has not been handled by AntTweakBar, process it
            {
                switch(event.type)
                {
                case sf::Event::Closed:
                    run = false;
                    break;
                default :
                    break;
                }
            }
        }
        //--------------------------------------------------------
        //--------------------Particle rendering------------------
        //--------------------------------------------------------

        //Bind shader
        glUseProgram(particleShader.id());
        ///Link uniform
        glUniformMatrix4fv(particleProjMatLoc, 1, GL_FALSE, glm::value_ptr(projMat));
        ///Link buffers and draw
        glEnableVertexAttribArray(particlePositionLoc);
        glBindBuffer(GL_ARRAY_BUFFER,particlePositionBuffer);
        glVertexAttribPointer(particlePositionLoc,2,GL_FLOAT,GL_FALSE,0,(void *)0);
        glEnableVertexAttribArray(particleColorLoc);
        glBindBuffer(GL_ARRAY_BUFFER,particleColorBuffer);
        glVertexAttribPointer(particleColorLoc,4,GL_FLOAT,GL_FALSE,0,(void *)0);
        glDrawArrays(GL_POINTS, 0, particle_center.size());
        glDisableVertexAttribArray(particlePositionLoc);
        glDisableVertexAttribArray(particleColorLoc);
        //Unbind shader
        glUseProgram(0);

        //--------------------------------------------------------
        //--------------------Grid rendering------------------
        //--------------------------------------------------------
        glUseProgram(gridShader.id());
        ///Link uniform
        glUniformMatrix4fv(gridProjMatLoc, 1, GL_FALSE, glm::value_ptr(projMat));
        ///Link buffers and draw
        glEnableVertexAttribArray(gridPositionLoc);
        glBindBuffer(GL_ARRAY_BUFFER,gridPositionBuffer);
        glVertexAttribPointer(gridPositionLoc,2,GL_FLOAT,GL_FALSE,0,(void *)0);
        glEnableVertexAttribArray(gridColorLoc);
        glBindBuffer(GL_ARRAY_BUFFER,gridColorBuffer);
        glVertexAttribPointer(gridColorLoc,4,GL_FLOAT,GL_FALSE,0,(void *)0);
        glDrawArrays(GL_LINE_LOOP, 0, 4);
        glDisableVertexAttribArray(gridPositionLoc);
        glDisableVertexAttribArray(gridColorLoc);
        glUseProgram(0);


        //--------------------------------------------------------
        //--------------------Fluid Surface rendering-------------
        //--------------------------------------------------------
        glUseProgram(fluidSurfaceShader.id());
        ///Link uniform
        glUniformMatrix4fv(fluidSurfaceProjMatLoc, 1, GL_FALSE, glm::value_ptr(projMat));
        ///Link buffers and draw
        glEnableVertexAttribArray(fluidSurfacePositionLoc);
        glBindBuffer(GL_ARRAY_BUFFER,fluidSurfacePositionBuffer);
        glVertexAttribPointer(fluidSurfacePositionLoc,2,GL_FLOAT,GL_FALSE,0,(void *)0);
        glEnableVertexAttribArray(fluidSurfaceColorLoc);
        glBindBuffer(GL_ARRAY_BUFFER,fluidSurfaceColorBuffer);
        glVertexAttribPointer(fluidSurfaceColorLoc,4,GL_FLOAT,GL_FALSE,0,(void *)0);
        glDrawArrays(GL_LINES, 0, 2.0*edges.size());
        glDisableVertexAttribArray(fluidSurfacePositionLoc);
        glDisableVertexAttribArray(fluidSurfaceColorLoc);
        glUseProgram(0);

        TwDraw();

        window.display();
    }

    //Deallocate buffers
    glDeleteBuffers(1, &particlePositionBuffer);
    glDeleteBuffers(1, &particleColorBuffer);
    glDeleteBuffers(1, &gridPositionBuffer);
    glDeleteBuffers(1, &fluidSurfacePositionBuffer);

    return 0;
}
