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
    int resolution = 1e3;
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
    sph.debugFluid();

    //Initialize visualization data
    std::vector< glm::vec2 > positions;
    std::vector< glm::vec4 > colors;
    std::vector< Vec2d > tmp_positions;

    tmp_positions = sph.getFluidPosition();
    for(size_t i=0; i<tmp_positions.size(); ++i)
    {
        positions.push_back(glm::vec2(tmp_positions[i][0], tmp_positions[i][1]));
        colors.push_back(glm::vec4(0.0,0.0,1.0,1.0));
    }
    tmp_positions = sph.getBoundaryPosition();
    for(size_t i=0; i<tmp_positions.size(); ++i)
    {
        positions.push_back(glm::vec2(tmp_positions[i][0], tmp_positions[i][1]));
        colors.push_back(glm::vec4(1.0,0.0,0.0,1.0));
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

    //Global settings
    glClearColor(0.0f,0.5f,0.5f,1.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_POINT_SPRITE);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    //Compile Shader and link shader
    Shader shader;
    string vertexPath = "./../particle_vertex.glsl";
    string fragmentPath = "./../particle_fragment.glsl";
    shader.load(vertexPath.c_str(), fragmentPath.c_str());

    //Set Shader variables
    GLuint positionLoc = glGetAttribLocation(shader.id(), "vertex_position");
    GLuint colorLoc = glGetAttribLocation(shader.id(), "vertex_color");

    //Allocate buffers
    GLuint positionBuffer;
    glGenBuffers(1, &positionBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, positionBuffer);
    glBufferData(GL_ARRAY_BUFFER, positions.size()*sizeof(glm::vec2), positions.data(), GL_STATIC_DRAW);

    GLuint colorBuffer;
    glGenBuffers(1, &colorBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, colorBuffer);
    glBufferData(GL_ARRAY_BUFFER, colors.size()*sizeof(glm::vec4), colors.data(), GL_STATIC_DRAW);


    bool run = true;
    while(run)
    {
        //Simulate
        sph.simulate();

        //Gather visualization data
        positions.clear();
        colors.clear();
        tmp_positions = sph.getFluidPosition();
        for(size_t i=0; i<tmp_positions.size(); ++i)
        {
            positions.push_back(glm::vec2(tmp_positions[i][0], tmp_positions[i][1]));
            colors.push_back(glm::vec4(0.0,0.0,1.0,1.0));
        }
        tmp_positions = sph.getBoundaryPosition();
        for(size_t i=0; i<tmp_positions.size(); ++i)
        {
            positions.push_back(glm::vec2(tmp_positions[i][0], tmp_positions[i][1]));
            colors.push_back(glm::vec4(1.0,0.0,0.0,1.0));
        }

        //Reload vbo
        glBindBuffer(GL_ARRAY_BUFFER, positionBuffer);
        glBufferData(GL_ARRAY_BUFFER, positions.size()*sizeof(glm::vec2), positions.data(), GL_STATIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, colorBuffer);
        glBufferData(GL_ARRAY_BUFFER, colors.size()*sizeof(glm::vec4), colors.data(), GL_STATIC_DRAW);

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

        //Bind shader
        glUseProgram(shader.id());

        ///Link buffers and draw
        glEnableVertexAttribArray(positionLoc);
        glBindBuffer(GL_ARRAY_BUFFER,positionBuffer);
        glVertexAttribPointer(positionLoc,2,GL_FLOAT,GL_FALSE,0,(void *)0);

        glEnableVertexAttribArray(colorLoc);
        glBindBuffer(GL_ARRAY_BUFFER,colorBuffer);
        glVertexAttribPointer(colorLoc,4,GL_FLOAT,GL_FALSE,0,(void *)0);

        glDrawArrays(GL_POINTS, 0, positions.size());

        glDisableVertexAttribArray(positionLoc);
        glDisableVertexAttribArray(colorLoc);

        //Unbind shader
        glUseProgram(0);

        window.display();
    }


    //Deallocate buffers
    glDeleteBuffers(1, &positionBuffer);
    glDeleteBuffers(1, &colorBuffer);

    return 0;
}
