#version 430

out vec4 colorOut; 

struct Materials {
	vec4 diffuse;
	vec4 ambient;
	vec4 specular;
	vec4 emissive;
	float shininess;
	int texCount;
};
uniform samplerCube cubeMap;
uniform bool normalMap;  //for normal mapping
uniform sampler2D normalMapForBumping;
uniform bool specularMap;
uniform uint diffMapCount;

uniform sampler2D fireworksTexMap;
uniform sampler2D treeTexMap;

uniform sampler2D texmap;
uniform sampler2D texmap1;
uniform sampler2D texmap2;
uniform sampler2D texMapFlare;
uniform sampler2D rollersTexMap;


uniform	sampler2D texUnitDiff;
uniform	sampler2D texUnitDiff1;
uniform	sampler2D texUnitSpec;
uniform	sampler2D texUnitNormalMap;

uniform int texMode;
uniform bool shadowMode;

uniform Materials mat;

uniform int headlights_switch;
uniform int pointlights_switch;
uniform int sunlight_switch;
uniform int fog_switch;

uniform vec4 dir_head;

vec4 diff,auxSpec;

#define NUM_TOTAL_LIGHTS 7
#define NUM_SPOTLIGHTS 2
#define POINT_L_ATTENUATION 0.3f
#define SPOT_L_ATTENUATION 0.7f
#define SUN_ATTENUATION 1.0f

in Data {
	vec3 normal;
	vec3 eye;
	vec3 lightDir;
	vec2 tex_coord;
} DataIn[NUM_TOTAL_LIGHTS];

in vec3 skyboxTexCoord;

const float reflect_factor = 0.9;


in vec4 pos;
uniform mat4 m_View;

void main() {

	vec4 spec = vec4(0.0);
	vec3 n, l, e;
	float intensity;

	vec4 accValue;
	vec4 texel = vec4(0.0);
	vec4 texel1 = vec4(0.0);
	vec4 cube_texel = vec4(0.0);
	colorOut = vec4(0.5, 0.5, 0.5, 1.0);

	bool a = false;
	 if (shadowMode)
	 	colorOut = vec4(0.5, 0.5, 0.5, 1.0); 
	 else {	
		for (int i = 0; i < NUM_TOTAL_LIGHTS-1; ++i) {

			// lookup normal from normal map, move from [0,1] to [-1, 1] range, normalize
			if(texMode == 6){  
				n = normalize(2.0 * texture(normalMapForBumping, DataIn[i].tex_coord).rgb - 1.0);
			} else {
				n = normalize(DataIn[i].normal);	
			}
			
			//If bump mapping, normalMap == TRUE, then lightDir and eye vectores come from vertex shader in tangent space
			vec3 l = normalize(DataIn[i].lightDir);
			vec3 e = normalize(DataIn[i].eye);

			float intensity = max(dot(n,l), 0.0);
			diff = mat.diffuse;
			auxSpec = mat.specular;
		
			if(i <= 3) {
				if (bool(pointlights_switch)) {
					vec3 l = normalize(DataIn[i].lightDir);
					vec3 e = normalize(DataIn[i].eye);

					float intensity = max(dot(n,l), 0.0);

					if (intensity > 0.0) {
						vec3 h = normalize(l + e);
						float intSpec = max(dot(h,n), 0.0);
						spec = auxSpec * pow(intSpec, mat.shininess);
						if(texMode == 1) {
							texel = texture(texmap, DataIn[i].tex_coord);
							texel1 = texture(texmap2, DataIn[i].tex_coord);
							accValue += (spec + intensity * mat.diffuse * texel * texel1) * POINT_L_ATTENUATION;
						} 
						else if(texMode == 2) {
							texel = texture(texUnitDiff, DataIn[i].tex_coord);
							accValue += (spec + intensity * diff * texel ) * POINT_L_ATTENUATION;
						} 
						else if (texMode == 4) { 
							// Billboard
							texel = texture(treeTexMap, DataIn[i].tex_coord);
							if (texel.a == 0.0) 
								discard;
							else {
								accValue += (spec + intensity * diff * texel) * POINT_L_ATTENUATION;
								accValue[3] += texel.a;
							}
						}
						else if (texMode == 5) { 
							//Fireworks
							texel = texture(fireworksTexMap, DataIn[i].tex_coord);
							if (texel.a == 0.0 || (mat.diffuse.a == 0.0)) 
								discard;
							else {
								accValue += (diff * texel) * POINT_L_ATTENUATION;
								accValue[3] += texel.a;
							}
						} 
						else if (texMode == 6) { 
							//rollers
							texel = texture(rollersTexMap, DataIn[i].tex_coord);
								accValue += (spec + intensity * mat.diffuse * texel ) * POINT_L_ATTENUATION;
						} 
						else if(texMode == 7) {
							// Environmental cube mapping
							vec3 reflected1 = vec3 (transpose(m_View) * vec4 (vec3(reflect(-e, n)), 0.0)); //reflection vector in world coord
							reflected1.x= -reflected1.x;   
							cube_texel = texture(cubeMap, reflected1);
							texel = texture(texmap, DataIn[i].tex_coord);  // texel from lighwood.tga
							vec4 aux_color = mix(texel, cube_texel, reflect_factor);
							aux_color = max(intensity*aux_color + spec, 0.1*aux_color);
							colorOut = vec4(aux_color.rgb, 1.0); 
						}
						else if(texMode == 9) {
								texel = texture(texmap, DataIn[i].tex_coord);
								texel1 = texture(texmap2, DataIn[i].tex_coord);
								accValue += vec4(max(intensity * mat.diffuse.rgb * texel.rgb * texel1.rgb + spec.xyz, mat.ambient.rgb * texel.rgb), 0.0);
						}
						else {
							accValue += (spec + intensity * diff) * POINT_L_ATTENUATION;
						}
					}
				}
			}
			else if (i <= 5) {
				if (bool(headlights_switch)) {

					float intensity = 0.0;
					vec4 spec = vec4(0.0);
					vec3 ld = normalize(DataIn[i].lightDir);

					float l_spotCutOff = 0.95f;
					vec4 l_spotDir = vec4(dir_head); // nao pode ser fixo

					vec3 sd = normalize(vec3(-l_spotDir));
					if (dot(sd,ld) > l_spotCutOff) {
						intensity = max(dot(n,ld), 0.0);
						if (intensity > 0.0) {
							vec3 eye = normalize(DataIn[i].eye);
							vec3 h = normalize(ld + eye);
							float intSpec = max(dot(h,n), 0.0);
							spec = auxSpec * pow(intSpec, mat.shininess);
							if(texMode == 1) {
								texel = texture(texmap, DataIn[i].tex_coord);
								texel1 = texture(texmap2, DataIn[i].tex_coord);
								accValue += (spec + intensity * mat.diffuse * texel * texel1) * SPOT_L_ATTENUATION;
							}
							else if(texMode == 2){
								texel = texture(texUnitDiff, DataIn[i].tex_coord);
								accValue += (spec + intensity * mat.diffuse * texel ) * SPOT_L_ATTENUATION;
							}
							else if (texMode == 4) { //Billboard
								texel = texture(treeTexMap, DataIn[i].tex_coord);
								if (texel.a == 0.0) 
									discard;
								else {
									accValue += (spec + intensity * diff * texel) * SPOT_L_ATTENUATION;
									accValue[3] += texel.a;
								}
							}
							else if (texMode == 5) { //Fireworks
								texel = texture(fireworksTexMap, DataIn[i].tex_coord);
								if (texel.a == 0.0 || (mat.diffuse.a == 0.0)) 
									discard;
								else {
									accValue += ( diff * texel) * SPOT_L_ATTENUATION;
									accValue[3] += texel.a;
								}
							}
							else if(texMode == 6) {
								texel = texture(rollersTexMap, DataIn[i].tex_coord);
								accValue += (spec + intensity * mat.diffuse * texel ) * SPOT_L_ATTENUATION;
							}
							else if(texMode == 7) {
								// Environmental cube mapping
								vec3 reflected1 = vec3 (transpose(m_View) * vec4 (vec3(reflect(-e, n)), 0.0)); //reflection vector in world coord
								reflected1.x= -reflected1.x;   
								cube_texel = texture(cubeMap, reflected1);
								texel = texture(texmap, DataIn[i].tex_coord);  // texel from lighwood.tga
								vec4 aux_color = mix(texel, cube_texel, reflect_factor);
								aux_color = max(intensity*aux_color + spec, 0.1*aux_color);
								colorOut = vec4(aux_color.rgb, 1.0); 
							}
							else if(texMode == 9) {
								texel = texture(texmap, DataIn[i].tex_coord);
								texel1 = texture(texmap2, DataIn[i].tex_coord);
								accValue += vec4(spec.rgb + intensity * mat.diffuse.rgb * texel.rgb * texel1.rgb, 0.0);
							}
							else {
								accValue += (spec + intensity * mat.diffuse) * SPOT_L_ATTENUATION;
							}
						}
					}
				}
			}
		} 

		if (bool(sunlight_switch)) {
			int dirLightIndex = NUM_TOTAL_LIGHTS - 1;
			vec4 spec = vec4(0.0);
	
			if(texMode == 6) {  // lookup normal from normal map, move from [0,1] to [-1, 1] range, normalize
				n = normalize(2.0 * texture(normalMapForBumping, DataIn[dirLightIndex].tex_coord).rgb - 1.0);
			}
			else{
				n = normalize(DataIn[dirLightIndex].normal);	
			}
			e = normalize(vec3(DataIn[dirLightIndex].eye));
			vec3 l_dir = vec3(0.0f, 0.0f, 1.0f);
	
			intensity = max(dot(n, l_dir), 0.0);
	
			if (intensity > 0.0) {
				vec3 h = normalize(l_dir + e);  
				float intSpec = max(dot(h,n), 0.0);
				spec = auxSpec * pow(intSpec, mat.shininess);
				if(texMode == 1) {
					texel = texture(texmap, DataIn[dirLightIndex].tex_coord);
					texel1 = texture(texmap2, DataIn[dirLightIndex].tex_coord);
					accValue += (spec + intensity * mat.diffuse * texel * texel1) * SUN_ATTENUATION;
				}
				else if(texMode == 2){
					texel = texture(texUnitDiff, DataIn[dirLightIndex].tex_coord);
					accValue += (spec + intensity * mat.diffuse * texel ) * SUN_ATTENUATION;
				}
				else if (texMode == 4) { //Billboard
					texel = texture(treeTexMap, DataIn[dirLightIndex].tex_coord);
					if (texel.a == 0.0) discard;
					else {
						accValue += (spec + intensity * diff * texel) * SUN_ATTENUATION;
						accValue[3] += texel.a;
					}
				}
				else if (texMode == 5) { //Fireworks
					texel = texture(fireworksTexMap, DataIn[dirLightIndex].tex_coord);
					if (texel.a == 0.0 || (mat.diffuse.a == 0.0)) discard;
					else {
						accValue += (diff * texel) * SUN_ATTENUATION;
						accValue[3] += texel.a;
					}
				}
				else if(texMode == 6){
					texel = texture(rollersTexMap, DataIn[dirLightIndex].tex_coord);
					accValue += (spec + intensity * mat.diffuse * texel ) * SUN_ATTENUATION;
				}
				else if(texMode == 7) // Environmental cube mapping
				{
					vec3 reflected1 = vec3 (transpose(m_View) * vec4 (vec3(reflect(-e, n)), 0.0)); //reflection vector in world coord
					reflected1.x= -reflected1.x;   
					cube_texel = texture(cubeMap, reflected1);
					texel = texture(texmap, DataIn[dirLightIndex].tex_coord);  // texel from lighwood.tga
					vec4 aux_color = mix(texel, cube_texel, reflect_factor);
					aux_color = max(intensity*aux_color + spec, 0.1*aux_color);
					accValue += vec4(aux_color.rgb, 1.0); 
				}

				else if(texMode == 9) {
								texel = texture(texmap, DataIn[dirLightIndex].tex_coord);
								texel1 = texture(texmap2, DataIn[dirLightIndex].tex_coord);
								accValue += vec4(spec.rgb + intensity * mat.diffuse.rgb * texel.rgb * texel1.rgb, 0.0);
							}

				else {
					accValue += (spec + intensity * mat.diffuse) * SUN_ATTENUATION;
				}

			}
		}		
		
		if(texMode == 1) { 
			colorOut = max(accValue, mat.ambient * texel * texel1);
		}
		else if (texMode == 2 || texMode == 4 || texMode == 5 || texMode == 6){
			colorOut = max(accValue, mat.ambient * texel);
		}
		else if(texMode == 3){
			vec4 texel = texture(cubeMap, skyboxTexCoord);
				colorOut = texel;
		}
		else if(texMode == 7){
			colorOut = accValue;
		}
		else if (texMode == 8) { //Lens Flare
			texel = texture(texMapFlare, DataIn[6].tex_coord);
			colorOut = texel * mat.diffuse;
		}
		else if(texMode == 9) {
			colorOut = max(accValue, mat.diffuse.a);
		}
		else {
			colorOut = max(accValue, mat.ambient);
		}

		// FOG ----------------------------------------------------------------------------

		if (bool(fog_switch)) {
			vec3 fogColor = vec3(0.255, 0.204, 0.204); 
			float d = length(pos); 
			float density = 0.05;
			float vis = exp(-d * density); 
			colorOut = mix(vec4(fogColor, 1.0), colorOut, vis);
		}


	 } // shadow mode


	// TRANSPARENCY
	// é definida no alpha channel dos arrays de propriedades de reflexao amb, diff, spec e emissive. 
	
}