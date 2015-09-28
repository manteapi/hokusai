#version 450

in vec4 fragmentColor;
out vec4 bufferColor;

void main()
{
    bufferColor = fragmentColor;
}
