#version 330 core

// Vertex attributes
layout(location = 0) in vec3 a_vertex;
layout(location = 1) in vec3 a_normal;

// Per-instance attributes
layout(location = 2) in vec4 i_color;
layout(location = 3) in mat4 i_matrix;

// Uniforms
uniform mat4 ModelViewMatrix;
uniform mat4 ProjectionMatrix;

// Outputs to the fragment shader
out vec3 v_position;
out vec3 v_normal;
out vec4 v_color;

void main() {
    mat4 NormalMatrix = inverse(transpose(ModelViewMatrix * i_matrix));
    vec4 worldPosition = i_matrix * vec4(a_vertex, 1.0);
    vec4 viewPosition  = ModelViewMatrix * worldPosition;    
    v_position = viewPosition.xyz;
    v_normal = (NormalMatrix * vec4(a_normal, 0)).xyz; 
    v_color = i_color;
    gl_Position = ProjectionMatrix * viewPosition;
}
