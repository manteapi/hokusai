#ifndef SHADER_HPP
#define SHADER_HPP

#include <GL/glew.h>
#include <string>

//debug info
void print_opengl_error(const char* file, int line);
void print_shader_info_log(GLuint shader);
void print_program_info_log(GLuint program);
std::string get_file_content(const std::string& filename);

class Shader
{
    public :
        Shader();
        ~Shader();
        void load(const char * vertex_file_path,const char * fragment_file_path);
        void reload(const char * vertex_file_path,const char * fragment_file_path);
        inline GLuint id() {return programId;}
    private :
        GLuint programId;
};

#endif
