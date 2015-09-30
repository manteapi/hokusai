#version 450

in vec2 fluidSurface_position;
in vec4 fluidSurface_color;
uniform mat4 projMat;
out vec4 fragmentColor;

void main()
{
    gl_PointSize=5.0f;
    gl_Position= projMat*vec4(fluidSurface_position, 0.0, 1.0);
    fragmentColor=fluidSurface_color;
}
