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

using namespace std;

int main()
{
    //Set simulation
    int resolution = 5e3;
    hokusai::System sph(resolution);

    ///Set fluid
    Vec2d fluidBox(0.5,0.5);
    Vec2d fluidOffset(0,0);
    sph.addParticleBox(fluidOffset, fluidBox);

    //Set boundary
    Vec2d securityOffset(1.05*sph.getSmoothingRadius());
    Vec2d boundBox(1.0,1.0);
    boundBox += securityOffset;
    Vec2d boundOffset(-0.5,-0.5);
    boundOffset -= 1.0*securityOffset;
    sph.addBoundaryBox(boundOffset, boundBox);

    ///Initialize parameters
    sph.init();
    sph.setTimeStep(1e-3);
    sph.debugFluid();
    sph.gridInfo.info();

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

    //Set Shader variables
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

        //Reload vbo                
        glBindBuffer(GL_ARRAY_BUFFER, particlePositionBuffer);
        glBufferData(GL_ARRAY_BUFFER, particle_center.size()*sizeof(glm::vec2), particle_center.data(), GL_STATIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, particleColorBuffer);
        glBufferData(GL_ARRAY_BUFFER, particle_color.size()*sizeof(glm::vec4), particle_color.data(), GL_STATIC_DRAW);

        ///Clear screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        ///Catch event
        sf::Event event;
        while(window.pollEvent(event))
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


        //--------------------------------------------------------
        //--------------------Particle rendering------------------
        //--------------------------------------------------------

        //Bind shader
        glUseProgram(particleShader.id());

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

        ///Link buffers and draw
        glEnableVertexAttribArray(gridPositionLoc);
        glBindBuffer(GL_ARRAY_BUFFER,gridPositionBuffer);
        glVertexAttribPointer(gridPositionLoc,2,GL_FLOAT,GL_FALSE,0,(void *)0);

        glEnableVertexAttribArray(gridColorLoc);
        glBindBuffer(GL_ARRAY_BUFFER,particleColorBuffer);
        glVertexAttribPointer(gridColorLoc,4,GL_FLOAT,GL_FALSE,0,(void *)0);

        glDrawArrays(GL_LINE_LOOP, 0, 4);

        glDisableVertexAttribArray(gridPositionLoc);
        glDisableVertexAttribArray(gridColorLoc);

        glUseProgram(0);

        window.display();
    }

    //Deallocate buffers
    glDeleteBuffers(1, &particlePositionBuffer);
    glDeleteBuffers(1, &particleColorBuffer);

    return 0;
}
