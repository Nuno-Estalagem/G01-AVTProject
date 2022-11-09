#version 430

uniform mat4 m_pvm;
uniform mat4 m_viewModel;
uniform mat3 m_normal;
uniform mat4 m_Model;   //por causa do cubo para a skybox


in vec4 position;
in vec4 normal;    //por causa do gerador de geometria
in vec4 texCoord;
uniform int texMode;
in vec4 tangent, bitangent;

#define NUM_POINT_LIGHTS 7

out Data {
	vec3 normal;
	vec3 eye;
	vec3 lightDir;
	vec2 tex_coord;

} DataOut[NUM_POINT_LIGHTS];

out vec3 skyboxTexCoord;

out vec4 pos;

uniform vec4 l_pos[NUM_POINT_LIGHTS]; // As posicoes das pointlights, definidas no lightDemo.cpp
uniform bool normalMap;

void main () {
	vec3 n, t, b;
	vec3 lightDir, eyeDir;
	vec3 aux;
	pos = m_viewModel * position;

	for (int i = 0; i < NUM_POINT_LIGHTS; ++i) {

		n = normalize(m_normal * normal.xyz);
		lightDir = vec3(l_pos[i] - pos);
		eyeDir = vec3(-pos);
		
		skyboxTexCoord = vec3(m_Model * position);	//Transformação de modelação do cubo unitário 
		skyboxTexCoord.x = - skyboxTexCoord.x; //Texturas mapeadas no interior logo negar a coordenada x
		DataOut[i].tex_coord = texCoord.st;

	if(texMode==6 || normalMap)  {  //transform eye and light vectors by tangent basis
			t = normalize(m_normal * tangent.xyz);
			b = tangent.w * cross(n,t);

			aux.x = dot(lightDir, t);
			aux.y = dot(lightDir, b);
			aux.z = dot(lightDir,n);
			lightDir = normalize(aux);

			aux.x = dot(eyeDir, t);
			aux.y = dot(eyeDir, b);
			aux.z = dot(eyeDir, n);
			eyeDir = normalize(aux);

		}
		DataOut[i].normal = n;
		DataOut[i].lightDir = lightDir;
		DataOut[i].eye = eyeDir;

	}

	

	gl_Position = m_pvm * position;	
}