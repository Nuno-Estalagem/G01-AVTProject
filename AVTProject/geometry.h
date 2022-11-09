#pragma once
#define MAX_TEXTURES 8

enum texType {DIFFUSE, SPECULAR, NORMALS, BUMP};

struct Material{
	float diffuse[4];
	float ambient[4];
	float specular[4];
	float emissive[4];
	float shininess;
	int texCount;
};
struct BetterMesh {
	float position[3] = { 0.0f };
	float size[3] = { 0.0f };
};
// A model can be made of many meshes. Each is stored  in the following structure
struct MyMesh : BetterMesh {
		GLuint vao;
		GLuint texUnits[MAX_TEXTURES];
		texType texTypes[4];
		float transform[16];
		GLuint numIndexes;
		unsigned int type;
		struct Material mat;
	};


MyMesh createCube();
MyMesh createCube(float size_x, float size_y, float size_z);
MyMesh createQuad(float size_x, float size_y);
MyMesh createSphere(float radius, int divisions);
MyMesh createTorus(float innerRadius, float outerRadius, int rings, int sides);
MyMesh createCylinder(float height, float radius, int sides);
MyMesh createCone(float height, float baseRadius, int sides);
MyMesh createPawn();
MyMesh computeVAO(int numP, float *p, float *pfloatoints, int sides, float smoothCos);
int revSmoothNormal2(float *p, float *nx, float *ny, float smoothCos, int beginEnd);
float *circularProfile(float minAngle, float maxAngle, float radius, int divisions, float transX= 0.0f, float transY = 0.0f);
void ComputeTangentArray(int vertexCount, float *vertex, float *normal, float *texcoord, GLuint indexesCount, GLuint *faceIndex, float *tangent);
