#version 450

in vec2 grid_position;
in vec4 grid_color;
out vec4 fragmentColor;

void main()
{
    gl_PointSize=5.0f;
    gl_Position= vec4(grid_position, 0.0, 1.0);
    fragmentColor=grid_color;
}
