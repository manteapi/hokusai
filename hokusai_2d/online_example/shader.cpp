#include <iostream>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include "shader.hpp"

using namespace std;

#define PRINT_OPENGL_ERROR() print_opengl_error(__FILE__, __LINE__);

void print_opengl_error(const char* file, int line)
{
    GLenum error;
    while( ( error=glGetError() ) != GL_NO_ERROR )
    {
        std::cerr << "glError in file" << file << ", line " << line << " : " << gluErrorString(error) << std::endl;
        error = glGetError();
    }

}

void print_shader_info_log(GLuint shader)
{
    int info_log_length = 0;
    int chars_written = 0;
    GLint status;

    glGetShaderiv(shader, GL_COMPILE_STATUS, &status); PRINT_OPENGL_ERROR();
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &info_log_length); PRINT_OPENGL_ERROR();

    if(info_log_length > 1)
    {
        std::vector<char> info_log(info_log_length+1);
        glGetShaderInfoLog(shader, info_log_length, &chars_written , &info_log[0]);
        std::cerr << "Shader InfoLog : " << std::endl << &info_log[0] << std::endl;
    }

    if( !status ) exit(EXIT_FAILURE);
}

void print_program_info_log(GLuint program)
{
    int     info_log_length = 0;
    int     chars_written  = 0;
    GLint   status;

    glGetProgramiv(program, GL_LINK_STATUS, &status); PRINT_OPENGL_ERROR();
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &info_log_length); PRINT_OPENGL_ERROR();

    if (info_log_length > 1)
    {
        std::vector<char> info_log(info_log_length+1);
        glGetProgramInfoLog(program, info_log_length, &chars_written, &info_log[0]); PRINT_OPENGL_ERROR ();
        std::cerr << "Program InfoLog : " << std::endl << &info_log[0] << std::endl;
    }
    if (!status) exit(EXIT_FAILURE);
}

std::string get_file_content(const std::string& filename)
{
    std::string code;
    std::ifstream stream(filename, std::ios::in);
    if(stream.is_open())
    {
        std::string Line = "";
        while(getline(stream, Line))
            code += "\n" + Line;
        stream.close();
    }
    else
    {
        std::cout << "Impossible d'ouvrir " << filename<< ". Are you in the right directory ?" << std::endl;
        exit(EXIT_FAILURE);
    }
    return code;
}

Shader::Shader() : programId(0) {}

Shader::~Shader()
{
    if(glIsProgram(programId))
    {
        glDeleteProgram(programId);
    }
}

void Shader::load(const char * vertex_file_path,const char * fragment_file_path)
{
    // Create & compile vertex shader 
    std::string VertexShaderCode = get_file_content(vertex_file_path);
    std::cout << "Compiling Shader : " << vertex_file_path << std::endl;
    char const * VertexSourcePointer = VertexShaderCode.c_str();
    GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL); PRINT_OPENGL_ERROR();
    glCompileShader(VertexShaderID); PRINT_OPENGL_ERROR();
    print_shader_info_log(VertexShaderID);

    // Create the fragment shader
    std::string FragmentShaderCode = get_file_content(fragment_file_path);
    std::cout << "Compiling shader : " << fragment_file_path << std::endl;
    char const * FragmentSourcePointer = FragmentShaderCode.c_str();
    GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL); PRINT_OPENGL_ERROR();
    glCompileShader(FragmentShaderID); PRINT_OPENGL_ERROR();
    print_shader_info_log(FragmentShaderID);

    //Create, attach, Link the program
    programId = glCreateProgram(); PRINT_OPENGL_ERROR();
    std::cout << "Linking program..." << std::endl;
    glAttachShader(programId, VertexShaderID); PRINT_OPENGL_ERROR();
    glAttachShader(programId, FragmentShaderID); PRINT_OPENGL_ERROR();
    glLinkProgram(programId); PRINT_OPENGL_ERROR();
    print_program_info_log(programId);

    //Delete vertex & fragment id
    glDeleteShader(VertexShaderID);
    glDeleteShader(FragmentShaderID);
}

void Shader::reload(const char * vertex_file_path,const char * fragment_file_path)
{
    //Check if program already contains a shader
    if(glIsProgram(programId))
    {
        glDeleteProgram(programId);
    }
    load(vertex_file_path, fragment_file_path);
}
