#version 450

in vec2 particle_position;
in vec4 particle_color;
uniform mat4 projMat;
out vec4 fragmentColor;

void main()
{
    gl_PointSize=3.0f;
    gl_Position= projMat*vec4(particle_position, 0.0, 1.0);
    fragmentColor=particle_color;
}
