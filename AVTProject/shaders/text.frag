#version 430

in vec2 TexCoords;
out vec4 color;

layout(binding = 0) uniform sampler2D text;
uniform vec3 textColor;

in vec3 position;
in vec3 normal, tangent, bitangent;
in vec2 texCoord;


void main()
{    
    vec4 sampled = vec4(1.0, 1.0, 1.0, texture(text, TexCoords).r);
    color = vec4(textColor, 1.0) * sampled;
}