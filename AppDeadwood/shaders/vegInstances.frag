#version 430 core

in vec3 v_position;
in vec3 v_normal;
in vec4 v_color;

out vec4 fragment;

// Phong shading components
uniform vec3 lightPos;
const   vec3 lightColor = vec3(1.0, 1.0, 1.0);

void main() {
    // Ambient component
    float ambientStrength = 0.2;
    vec3 ambient = ambientStrength * lightColor;

    // Diffuse component
    vec3 norm = normalize(v_normal);
    vec3 lightDir = normalize(lightPos - v_position);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = (1 - ambientStrength) * diff * lightColor;

    // Specular component
    float specularStrength = 0.5;
    vec3 viewDir = normalize(-v_position);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32.0);
    vec3 specular = specularStrength * spec * lightColor;

    // Combine components
    vec3 result = (ambient + diffuse)*v_color.rgb + specular;
    fragment = vec4(result, v_color.a);
}
