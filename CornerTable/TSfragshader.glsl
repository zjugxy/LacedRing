#version 450 core
layout(location = 0) out vec4 FragColor;

in myPerVertexData
{
    vec3 color;
} fragIn; 

void main()
{
  FragColor = vec4(fragIn.color,1.0);
}