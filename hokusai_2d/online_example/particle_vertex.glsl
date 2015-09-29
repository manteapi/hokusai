#version 450

in vec2 particle_position;
in vec4 particle_color;
out vec4 fragmentColor;

void main()
{
    gl_PointSize=10.0f;
    gl_Position= vec4(particle_position, 0.0, 1.0);
    fragmentColor=particle_color;
}
