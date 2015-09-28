#version 450

in vec2 vertex_position;
in vec4 vertex_color;

out vec4 fragmentColor;

void main()
{
    gl_PointSize=5.0f;
    gl_Position= vec4(vertex_position, 0.0, 1.0);
    fragmentColor=vertex_color;
}
