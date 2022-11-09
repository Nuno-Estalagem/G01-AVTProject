//
// AVT: Phong Shading and Text rendered with FreeType library
// The text rendering was based on https://learnopengl.com/In-Practice/Text-Rendering
// This demo was built for learning purposes only.
// Some code could be severely optimised, but I tried to
// keep as simple and clear as possible.
//
// The code comes with no warranties, use it at your own risk.
// You may use it, or parts of it, wherever you want.
//
// Author: Jo�o Madeiras Pereira
//

#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>

// include GLEW to access OpenGL 3.3 functions
#include <GL/glew.h>

// GLUT is the toolkit to interface with the OS
#include <GL/freeglut.h>

#include <IL/il.h>

#include "assimp/Importer.hpp" //OO version Header!
#include "assimp/scene.h"

// Use Very Simple Libs
#include "VSShaderlib.h"
#include "AVTmathLib.h"
#include "VertexAttrDef.h"
#include "geometry.h"
#include "Rover.h"
#include "RandomRoller.h"
#include "avtFreeType.h"
#include "meshFromAssimp.h"

#include "Camera.h"
#include <Texture_Loader.h>

#include "l3dBillboard.h"
#include <flare.h>

#define frand() ((float)rand() / RAND_MAX)
#define M_PI 3.14159265
#define MAX_PARTICULAS 1500

using namespace std;

#define CAPTION "AVT Demo: Phong Shading and Text rendered with FreeType"
int WindowHandle = 0;
int WinX = 1024, WinY = 768;

unsigned int FrameCount = 0;

// Created an instance of the Importer class in the meshFromAssimp.cpp file
extern Assimp::Importer importer;
// the global Assimp scene object
extern const aiScene *scene;

// shaders
VSShaderLib shader;		// geometry
VSShaderLib shaderText; // render bitmap text

// File with the font
const string font_name = "fonts/arial.ttf";

// Vector with meshes
vector<struct MyMesh> spyMeshes;
vector<struct MyMesh> myMeshes;
std::vector<struct MyMesh> billboardMesh;
std::vector<struct MyMesh> skyBoxMesh;
std::vector<struct MyMesh> recreativeMeshes;
std::vector<struct MyMesh> flareMeshes;

// External array storage defined in AVTmathLib.cpp

/// The storage for matrices
extern float mMatrix[COUNT_MATRICES][16];
extern float mCompMatrix[COUNT_COMPUTED_MATRICES][16];

/// The normal matrix
extern float mNormal3x3[9];

GLint pvm_uniformId;
GLint vm_uniformId;
GLint normal_uniformId;
GLint lPos_uniformId;
GLint tex_loc, tex_loc1, tex_loc2, tex_loc_cube, tex_sphereMap_loc, tex_normalMap_loc, view_uniformId;
GLint tex_loc3, tex_loc4, tex_loc5, tex_loc6;

GLint lPos_uniformId_1;
GLint lPos_uniformId_2;
GLint lPos_uniformId_3;
GLint lPos_uniformId_4;
GLint lPos_uniformId_5;
GLint lPos_uniformId_6;
GLint lPos_uniformId_7;
GLint lPos_uniformId_8;
GLint lPos_uniformId_9;

GLint model_uniformId;
GLint shadowMode_uniformId;
GLint dir_headlight_uniformId;
GLint headlights_uniformId;
GLint pointlights_uniformId;
GLint sunlight_uniformId;
GLint fog_uniformId;
GLint texMode_uniformId;
GLuint FlareTextureArray[5];
GLuint TextureArray[23];
bool normalMapKey = TRUE;

GLint normalMap_loc;
GLint specularMap_loc;
GLint diffMapCount_loc;
GLint reflect_perFragment_uniformId;

float billboard_type = 0;
float cam_X, cam_Y, cam_Z;

typedef struct
{
	float life;			// vida
	float fade;			// fade
	float r, g, b;		// color
	GLfloat x, y, z;	// posi��o
	GLfloat vx, vy, vz; // velocidade
	GLfloat ax, ay, az; // acelera��o
} Particle;

Particle particula[MAX_PARTICULAS];
int dead_num_particles = 0;
int fireworks = 0;

int objId = 0;
// std::vector<struct MyMesh> mesh;

// Camera Position
// float camX, camY, camZ;

// Mouse Tracking Variables
int startX, startY, tracking = 0;

// Camera Spherical Coordinates
float alpha = 39.0f, beta = 51.0f;
float r = 10.0f;

// Frame counting and FPS computation
long myTime, timebase = 0, frame = 0;
char s[32];
// float lightPos[4] = { 4.0f, 6.0f, 2.0f, 1.0f };

float floorX = 200.0f;
float floorY = 0.005f;
float floorZ = 200.0f;

float headlight_1_pos[4];
float headlight_2_pos[4];

Rover rover;
const int numRollers = 4;
RandomRoller randomRollers[numRollers];
const float rangeX = floorX;
const float rangeZ = floorZ;
bool moveRollers = true;

// acceleration
float acceleration_x = 1.0f;
float velocity_x = 0.0f;
float position_x = 0.0f;
float acceleration_z = 15.0f;
float velocity_z = 0.0f;
float position_z = 0.0f;
float deltaT = 1.0f / 50.0f;
const float acceleration_down_x = 10;
const float acceleration_up_x = 15;
const float acceleration_down_z = 25;
const float acceleration_up_z = 30;
const float velocity_max_x = 20;
const float velocity_max_z = 30;
const float velocity_min = 0.01f;

enum MovementState
{
	IDLE,
	FORWARD,
	BACKWARD,
	LEFT,
	RIGHT
};
char *MovementType[] = {
	"IDLE",
	"FORWARD",
	"BACKWARD",
	"LEFT",
	"RIGHT"};
enum PauseState
{
	PAUSE,
	PLAY,
	TRANSITION_PAUSE,
	TRANSITION_PLAY,
	GAME_FINISHED,
	TRANSITION_FINISHED
};
char *PauseType[] = {
	"PAUSE",
	"PLAY",
	"TRANSITION_PAUSE",
	"TRANSITION_PLAY"
	"GAME_FINISHED"};
PauseState gameState = PLAY;
MovementState movStateX = IDLE;
MovementState movStateZ = IDLE;

int points = 0;
int baseNumLives = 1;
int numLives = baseNumLives;
char model_dir[50] = {""};

bool movementKey[] = {
	false, // forward
	false, // right
	false, // backward
	false  // left
};

bool wasUpx = false;
bool wasUpz = false;

ofstream test_file;

// NEW CAMERA CODE ------------------

int active_camera = 0;
Camera cams[3];

float ratio = 1024.0f / 768.0f;
float cam_dist = -5.0f;
float cam_height = 3.0f;
float cam_offset;

const int NUM_POINT_LIGHTS = 7;

float pointLightPositions[NUM_POINT_LIGHTS][4] = {
	{15.0f, 10.0f, 0.0f, 1.0f},	  // pointlight (1 of 6)
	{0.0f, 10.0f, 15.0f, 1.0f},	  // pointlight (2 of 6)
	{-15.0f, 10.0f, 0.0f, 1.0f},  // pointlight (3 of 6)
	{0.0f, 10.0f, -15.0f, 1.0f},  // pointlight (4 of 6)
	{-15.0f, 2.0f, -15.0f, 1.0f}, // pointlight (5 of 6)
	{15.0f, 9.0f, 15.0f, 1.0f},	  // pointlight (6 of 6)
	{0.0f, -1.0f, 0.0f, 0.0f}	  // directional light (x, y (altura), z)
};

bool headlights = true;
bool pointlights = true;
bool sunlight = true;
bool fog = true;

float lightScreenPos[3];
FLARE_DEF AVTflare;

bool changedSize = false;

int reflect_perFragment = 0;

float pause_velocity_x;
float pause_velocity_z;

MyMesh rearCubeCamera;

bool shadowMode = false;

inline double clamp(const double x, const double min, const double max)
{
	return (x < min ? min : (x > max ? max : x));
}

inline int clampi(const int x, const int min, const int max)
{
	return (x < min ? min : (x > max ? max : x));
}

unsigned int getTextureId(char *name)
{
	int i;

	for (i = 0; i < NTEXTURES; ++i)
	{
		if (strncmp(name, flareTextureNames[i], strlen(name)) == 0)
			return i;
	}
	return -1;
}

void loadFlareFile(FLARE_DEF *flare, char *filename)
{
	int n = 0;
	FILE *f;
	char buf[256];
	int fields;

	memset(flare, 0, sizeof(FLARE_DEF));

	f = fopen(filename, "r");
	if (f)
	{
		fgets(buf, sizeof(buf), f);
		sscanf(buf, "%f %f", &flare->fScale, &flare->fMaxSize);

		while (!feof(f))
		{
			char name[8] = {
				'\0',
			};
			double dDist = 0.0, dSize = 0.0;
			float color[4];
			int id;

			fgets(buf, sizeof(buf), f);
			fields = sscanf(buf, "%4s %lf %lf ( %f %f %f %f )", name, &dDist, &dSize, &color[3], &color[0], &color[1], &color[2]);
			if (fields == 7)
			{
				for (int i = 0; i < 4; ++i)
					color[i] = clamp(color[i] / 255.0f, 0.0f, 1.0f);
				id = getTextureId(name);
				if (id < 0)
					printf("Texture name not recognized\n");
				else
					flare->element[n].textureId = id;
				flare->element[n].fDistance = (float)dDist;
				flare->element[n].fSize = (float)dSize;
				memcpy(flare->element[n].matDiffuse, color, 4 * sizeof(float));
				++n;
			}
		}

		flare->nPieces = n;
		fclose(f);
	}
	else
		printf("Flare file opening error\n");
}

void iniParticles(void)
{
	GLfloat v, theta, phi;
	int i;

	for (i = 0; i < MAX_PARTICULAS; i++)
	{
		v = 0.8 * frand() + 0.2;
		phi = frand() * M_PI;
		theta = 2.0 * frand() * M_PI;

		particula[i].x = 0.0f;
		particula[i].y = 10.0f;
		particula[i].z = 0.0f;
		particula[i].vx = v * cos(theta) * sin(phi);
		particula[i].vy = v * cos(phi);
		particula[i].vz = v * sin(theta) * sin(phi);
		particula[i].ax = 0.1f;	  /* simular um pouco de vento */
		particula[i].ay = -0.15f; /* simular a acelera��o da gravidade */
		particula[i].az = 0.0f;

		/* tom amarelado que vai ser multiplicado pela textura que varia entre branco e preto */
		particula[i].r = 0.882f;
		particula[i].g = 0.552f;
		particula[i].b = 0.211f;

		particula[i].life = 1.0f;	 /* vida inicial */
		particula[i].fade = 0.0025f; /* step de decr�scimo da vida para cada itera��o */
	}
}

void updateParticles()
{
	int i;
	float h;

	/* M�todo de Euler de integra��o de eq. diferenciais ordin�rias
	h representa o step de tempo; dv/dt = a; dx/dt = v; e conhecem-se os valores iniciais de x e v */

	// h = 0.125f;
	h = 0.033;
	if (fireworks && gameState!=PAUSE)
	{

		for (i = 0; i < MAX_PARTICULAS; i++)
		{
			particula[i].x += (h * particula[i].vx);
			particula[i].y += (h * particula[i].vy);
			particula[i].z += (h * particula[i].vz);
			particula[i].vx += (h * particula[i].ax);
			particula[i].vy += (h * particula[i].ay);
			particula[i].vz += (h * particula[i].az);
			particula[i].life -= particula[i].fade;
		}
	}
}

void drawTreeBillboards()
{
	// Draw trees billboards. Despite of blending, they are considered opaque so z-buffer works normally
	GLint loc;
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	float pos[3];
	float cam[3] = {cam_X, cam_Y, cam_Z};

	// glUniform1i(texMode_uniformId, 1); // draw textured quads
	glUniform1i(texMode_uniformId, 4);

	for (int i = -5; i < 5; i++)
		for (int j = -5; j < 5; j++)
		{
			pushMatrix(MODEL);
			translate(MODEL, 5 + i * 10.0, 0, 5 + j * 10.0);

			pos[0] = 5 + i * 10.0;
			pos[1] = 0;
			pos[2] = 5 + j * 10.0;

			/*
			cout << "pos[0]" << pos[0] << endl;
			cout << "pos[1]" << pos[1] << endl;
			cout << "pos[2]" << pos[2] << endl;
			cout << "--------------------" << endl;*/

			if (billboard_type == 2)
				l3dBillboardSphericalBegin(cam, pos);
			else if (billboard_type == 3)
				l3dBillboardCylindricalBegin(cam, pos);

			objId = 0; // quad for tree

			// diffuse and ambient color are not used in the tree quads
			loc = glGetUniformLocation(shader.getProgramIndex(), "mat.ambient");
			glUniform4fv(loc, 1, billboardMesh[objId].mat.ambient);
			loc = glGetUniformLocation(shader.getProgramIndex(), "mat.diffuse");
			glUniform4fv(loc, 1, billboardMesh[objId].mat.diffuse);
			loc = glGetUniformLocation(shader.getProgramIndex(), "mat.specular");
			glUniform4fv(loc, 1, billboardMesh[objId].mat.specular);
			loc = glGetUniformLocation(shader.getProgramIndex(), "mat.shininess");
			glUniform1f(loc, billboardMesh[objId].mat.shininess);

			pushMatrix(MODEL);
			translate(MODEL, 0.0, 3.0, 0.0f);

			// send matrices to OGL
			if (billboard_type == 0 || billboard_type == 1)
			{ // Cheating matrix reset billboard techniques
				computeDerivedMatrix(VIEW_MODEL);

				// reset VIEW_MODEL
				if (billboard_type == 0)
					BillboardCheatSphericalBegin();
				else
					BillboardCheatCylindricalBegin();

				computeDerivedMatrix_PVM(); // calculate PROJ_VIEW_MODEL
			}
			else
				computeDerivedMatrix(PROJ_VIEW_MODEL);

			glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
			glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
			computeNormalMatrix3x3();
			glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);
			glBindVertexArray(billboardMesh[objId].vao);
			glDrawElements(billboardMesh[objId].type, billboardMesh[objId].numIndexes, GL_UNSIGNED_INT, 0);
			popMatrix(MODEL);

			//	if (type==0 || type==1) // restore matrix VIEW_MODEL n�o � necess�rio pois a PVM � sempre calculada a pArtir da MODEL e da VIEW que n�o s�o ALTERADAS

			popMatrix(MODEL);
		}

	if (fireworks)
	{

		float particle_color[4];

		updateParticles();

		// draw fireworks particles
		objId = 1; // quad for particle

		glUniform1i(texMode_uniformId, 5); // particle.tga

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glDepthMask(GL_FALSE); // Depth Buffer Read Only

		// glUniform1i(texMode_uniformId, 2); // draw modulated textured particles
		//  ISTO PODE ESTAR BUGADEX, PORQUE NAO ATUALIZEI ISTO -------------------------------------------------------------------

		for (int i = 0; i < MAX_PARTICULAS; i++)
		{
			if (particula[i].life > 0.0f) /* s� desenha as que ainda est�o vivas */
			{

				// cout << "this is a live particle" << endl;

				/* A vida da part�cula representa o canal alpha da cor. Como o blend est� activo a cor final � a soma da cor rgb do fragmento multiplicada pelo
				alpha com a cor do pixel destino */

				particle_color[0] = particula[i].r;
				particle_color[1] = particula[i].g;
				particle_color[2] = particula[i].b;
				particle_color[3] = particula[i].life;

				// send the material - diffuse color modulated with texture
				loc = glGetUniformLocation(shader.getProgramIndex(), "mat.diffuse");
				glUniform4fv(loc, 1, particle_color);

				pushMatrix(MODEL);
				translate(MODEL, particula[i].x, particula[i].y, particula[i].z);

				// send matrices to OGL
				computeDerivedMatrix(PROJ_VIEW_MODEL);
				glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
				glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
				computeNormalMatrix3x3();
				glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

				glBindVertexArray(billboardMesh[objId].vao);
				glDrawElements(billboardMesh[objId].type, billboardMesh[objId].numIndexes, GL_UNSIGNED_INT, 0);
				popMatrix(MODEL);
			}
			else
				dead_num_particles++;
		}

		glDepthMask(GL_TRUE); // make depth buffer again writeable

		if (dead_num_particles == MAX_PARTICULAS)
		{
			fireworks = 0;
			dead_num_particles = 0;
			printf("All particles dead\n");
		}
	}
	glBindVertexArray(0);
	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_BLEND);
}

void render_flare(FLARE_DEF *flare, int lx, int ly, int *m_viewport)
{ // lx, ly represent the projected position of light on viewport

	int dx, dy; // Screen coordinates of "destination"
	int px, py; // Screen coordinates of flare element
	int cx, cy;
	float maxflaredist, flaredist, flaremaxsize, flarescale, scaleDistance;
	int width, height, alpha; // Piece parameters;
	int i;
	float diffuse[4];

	GLint loc;

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	int screenMaxCoordX = m_viewport[0] + m_viewport[2] - 1;
	int screenMaxCoordY = m_viewport[1] + m_viewport[3] - 1;

	// viewport center
	cx = m_viewport[0] + (int)(0.5f * (float)m_viewport[2]) - 1;
	cy = m_viewport[1] + (int)(0.5f * (float)m_viewport[3]) - 1;

	// Compute how far off-center the flare source is.
	maxflaredist = sqrt(cx * cx + cy * cy);
	flaredist = sqrt((lx - cx) * (lx - cx) + (ly - cy) * (ly - cy));
	scaleDistance = (maxflaredist - flaredist) / maxflaredist;
	flaremaxsize = (int)(m_viewport[2] * flare->fMaxSize);
	flarescale = (int)(m_viewport[2] * flare->fScale);

	// Destination is opposite side of centre from source
	dx = clampi(cx + (cx - lx), m_viewport[0], screenMaxCoordX);
	dy = clampi(cy + (cy - ly), m_viewport[1], screenMaxCoordY);

	// Render each element. To be used Texture Unit 0

	glUniform1i(texMode_uniformId, 8); // draw modulated textured particles
	glUniform1i(tex_loc6, 0);		   // use TU 0

	for (i = 0; i < flare->nPieces; ++i)
	{
		// Position is interpolated along line between start and destination.
		px = (int)((1.0f - flare->element[i].fDistance) * lx + flare->element[i].fDistance * dx);
		py = (int)((1.0f - flare->element[i].fDistance) * ly + flare->element[i].fDistance * dy);
		px = clampi(px, m_viewport[0], screenMaxCoordX);
		py = clampi(py, m_viewport[1], screenMaxCoordY);

		// Piece size are 0 to 1; flare size is proportion of screen width; scale by flaredist/maxflaredist.
		width = (int)(scaleDistance * flarescale * flare->element[i].fSize);

		// Width gets clamped, to allows the off-axis flaresto keep a good size without letting the elements get big when centered.
		if (width > flaremaxsize)
			width = flaremaxsize;

		height = (int)((float)m_viewport[3] / (float)m_viewport[2] * (float)width);
		memcpy(diffuse, flare->element[i].matDiffuse, 4 * sizeof(float));
		diffuse[3] *= scaleDistance; // scale the alpha channel

		if (width > 1)
		{
			// send the material - diffuse color modulated with texture
			loc = glGetUniformLocation(shader.getProgramIndex(), "mat.diffuse");
			glUniform4fv(loc, 1, diffuse);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, FlareTextureArray[flare->element[i].textureId]);
			pushMatrix(MODEL);
			translate(MODEL, (float)(px - width * 0.0f), (float)(py - height * 0.0f), 0.0f);
			scale(MODEL, (float)width, (float)height, 1);
			computeDerivedMatrix(PROJ_VIEW_MODEL);
			glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
			glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
			computeNormalMatrix3x3();
			glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

			glBindVertexArray(flareMeshes[0].vao);
			glDrawElements(flareMeshes[0].type, flareMeshes[0].numIndexes, GL_UNSIGNED_INT, 0);
			glBindVertexArray(0);
			popMatrix(MODEL);
		}
	}
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glDisable(GL_BLEND);
}

void showFlares()
{

	if (sunlight && active_camera == 2 && !fog)
	{

		int flarePos[2];
		int m_viewport[4];
		glGetIntegerv(GL_VIEWPORT, m_viewport);

		pushMatrix(MODEL);
		loadIdentity(MODEL);
		computeDerivedMatrix(PROJ_VIEW_MODEL); // pvm to be applied to lightPost. pvm is used in project function
		if (!project(pointLightPositions[3], lightScreenPos, m_viewport))
			printf("Error in getting projected light in screen\n"); // Calculate the window Coordinates of the light position: the projected position of light on viewport
		flarePos[0] = clampi((int)lightScreenPos[0], m_viewport[0], m_viewport[0] + m_viewport[2] - 1);
		flarePos[1] = clampi((int)lightScreenPos[1], m_viewport[1], m_viewport[1] + m_viewport[3] - 1);
		popMatrix(MODEL);

		// viewer looking down at  negative z direction
		pushMatrix(PROJECTION);
		loadIdentity(PROJECTION);
		pushMatrix(VIEW);
		loadIdentity(VIEW);
		ortho(m_viewport[0], m_viewport[0] + m_viewport[2] - 1, m_viewport[1], m_viewport[1] + m_viewport[3] - 1, -1, 1);
		render_flare(&AVTflare, flarePos[0], flarePos[1], m_viewport);
		popMatrix(PROJECTION);
		popMatrix(VIEW);
	}
}

float randf()
{
	return ((float)rand()) / ((float)RAND_MAX);
}

void timer(int value)
{
	std::ostringstream oss;
	oss << CAPTION << ": " << FrameCount << " FPS @ (" << WinX << "x" << WinY << ")";
	std::string s = oss.str();
	glutSetWindow(WindowHandle);
	glutSetWindowTitle(s.c_str());
	FrameCount = 0;
	glutTimerFunc(1000, timer, 0);
}

void refresh(int value)
{
	glutPostRedisplay();
	glutTimerFunc(1000 / value, refresh, value);
}

// ------------------------------------------------------------
//
// Reshape Callback Function
//

void changeSize(int w, int h)
{
	changedSize = true;
	float ratio;
	// Prevent a divide by zero, when window is too short
	if (h == 0)
		h = 1;
	// set the viewport to be the entire window
	glViewport(0, 0, w, h);

	/* create a diamond shaped stencil area */
	loadIdentity(PROJECTION);

	//// set the projection matrix
	ratio = (1.0f * w) / h;
	loadIdentity(PROJECTION);
	perspective(53.13f, ratio, 0.1f, 1000.0f);

	WinX = w;
	WinY = h;

	// ortho(0.0f, 1024, 0.0f, 768, -10, 10
	// ortho(-32.0 * (GLfloat)w / (GLfloat)h, 32.0 * (GLfloat)w / (GLfloat)h,-64.0, 64.0, -10, 10
}

// ------------------------------------------------------------
//
// Render stuff

// It setups the camera into the desired mode, which in the last one
// makes the camera following the rover
void setupCamera(bool beforeStencil)
{
	const float *roverCenter = rover.getRoverCenter();

	if (active_camera == 0)
	{
		ortho(-floorX / 2, floorX / 2, -floorZ / 2, floorZ / 2, -1000, 1000);

		lookAt(cams[active_camera].camPos[0] + 0.01f, cams[active_camera].camPos[1], cams[active_camera].camPos[2],
			   cams[active_camera].camTarget[0], cams[active_camera].camTarget[1], cams[active_camera].camTarget[2],
			   0, 1, 0);
	}
	else if (active_camera == 1)
	{
		perspective(150.13f, ratio, 0.1f, 1000.0f);
		lookAt(cams[active_camera].camPos[0] + 0.01f, cams[active_camera].camPos[1], cams[active_camera].camPos[2],
			   cams[active_camera].camTarget[0], cams[active_camera].camTarget[1], cams[active_camera].camTarget[2],
			   0, 1, 0);
	}
	else if (active_camera == 2)
	{
		// Before
		float camaraPositionX = cams[active_camera].camPos[0] + roverCenter[0];
		float camaraPositionZ = cams[active_camera].camPos[2] + roverCenter[2];
		float camaraPositionY = cams[active_camera].camPos[1] + roverCenter[1];

		cams[active_camera].camTarget[0] = roverCenter[0];
		cams[active_camera].camTarget[1] = roverCenter[1];
		cams[active_camera].camTarget[2] = roverCenter[2];

		perspective(53.13f, ratio, 0.1f, 10000.0f);

		lookAt(camaraPositionX + 0.01f, camaraPositionY, camaraPositionZ,
				cams[active_camera].camTarget[0], cams[active_camera].camTarget[1], cams[active_camera].camTarget[2],
				0, 1, 0);

		translate(VIEW, roverCenter[0], roverCenter[1], roverCenter[2]);
		rotate(VIEW, -*rover.getRoverAngle(), 0.0, 1.0, 0.0);
		translate(VIEW, -roverCenter[0], -roverCenter[1], -roverCenter[2]);
	}
	else if (active_camera == 3) {
		// todo
		const float* rearCameraConfig = rover.getRearCameraConfiguration();
		if (beforeStencil) {

			//perspective(83.13f, (1.0f * WinX) / WinY, 0.1f, 1000.0f);
			perspective(63.13f, ratio, 0.1f, 1000.0f);
			//scale(VIEW, -1, 1, 1);
			lookAt(
				rearCameraConfig[0],
				rearCameraConfig[1],
				rearCameraConfig[2],
				rearCameraConfig[3],
				rearCameraConfig[4],
				rearCameraConfig[5],
				0,
				1,
				0);
		}
		else {
			perspective(53.13f, ratio, 0.1f, 1000.0f);
			lookAt(
				rearCameraConfig[6],
				rearCameraConfig[7],
				rearCameraConfig[8],
				rearCameraConfig[9],
				rearCameraConfig[10],
				rearCameraConfig[11],
				0,
				1,
				0);
		}
	}
}

// glMultMatrix aux function aType = aMatrix * aType
void multMatrixWithType(float* aMatrix, MatrixTypes aType)
{

	float* a, * b, res[16];
	a = aMatrix;
	b = mMatrix[aType];

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			res[j * 4 + i] = 0.0f;
			for (int k = 0; k < 4; ++k) {
				res[j * 4 + i] += a[k * 4 + i] * b[j * 4 + k];
			}
		}
	}
	memcpy(mMatrix[aType], res, 16 * sizeof(float));
}


// Function to perform the movement on the object in the scene
void moveSceneObjects()
{
	for (int i = 0; i < numRollers; i++)
	{
		randomRollers[i].drawRoller();
	}
	rover.renderRover();
}

// Setup the directional, spotlights and head lights
void setupLightInScene()
{
	// send the light position in eye coordinates
	// glUniform4fv(lPos_uniformId, 1, lightPos); //efeito capacete do mineiro, ou seja lighPos foi definido em eye coord

	float res[4];
	glUniform1i(headlights_uniformId, headlights);
	glUniform1i(pointlights_uniformId, pointlights);
	glUniform1i(sunlight_uniformId, sunlight);
	glUniform1i(fog_uniformId, fog);

	float res1[4];
	multMatrixPoint(VIEW, pointLightPositions[0], res1); // NEW
	glUniform4fv(lPos_uniformId_1, 1, res1);

	float res2[4];
	multMatrixPoint(VIEW, pointLightPositions[1], res2); // NEW
	glUniform4fv(lPos_uniformId_2, 1, res2);

	float res5[4];
	multMatrixPoint(VIEW, pointLightPositions[2], res5); // NEW
	glUniform4fv(lPos_uniformId_3, 1, res5);

	float res6[4];
	multMatrixPoint(VIEW, pointLightPositions[3], res6); // NEW
	glUniform4fv(lPos_uniformId_4, 1, res6);

	float res7[4];
	multMatrixPoint(VIEW, pointLightPositions[4], res7); // NEW
	glUniform4fv(lPos_uniformId_5, 1, res7);

	float res8[4];
	multMatrixPoint(VIEW, pointLightPositions[5], res8); // NEW
	glUniform4fv(lPos_uniformId_6, 1, res8);

	float *headLighterPos = (float *)rover.getHeadLightConfiguration();

	headlight_1_pos[0] = headLighterPos[0]; headlight_2_pos[0] = headLighterPos[3];
	headlight_1_pos[1] = headLighterPos[1]; headlight_2_pos[1] = headLighterPos[4];
	headlight_1_pos[2] = headLighterPos[2]; headlight_2_pos[2] = headLighterPos[5];

	float headLightDirection[3] = { headLighterPos[6],headLighterPos[7],headLighterPos[8] };

	float res3[4];
	multMatrixPoint(VIEW, headlight_1_pos, res3); // SP2
	glUniform4fv(lPos_uniformId_7, 1, res3);

	float res4[4];
	multMatrixPoint(VIEW, headlight_2_pos, res4); // SP1
	glUniform4fv(lPos_uniformId_8, 1, res4);

	float res10[4];
	multMatrixPoint(VIEW, headLightDirection, res10);
	glUniform4fv(dir_headlight_uniformId, 1, res10);

	//float res9[4];
	//multMatrixPoint(VIEW, pointLightPositions[6], res9); // DIRECTIONAL LIGHT
	//glUniform4fv(lPos_uniformId_9, 1, res9);
}

// Setupping the textures
void setupTextures()
{
	glActiveTexture(GL_TEXTURE4);
	glBindTexture(GL_TEXTURE_2D, TextureArray[0]);

	glActiveTexture(GL_TEXTURE5);
	glBindTexture(GL_TEXTURE_2D, TextureArray[1]);

	glActiveTexture(GL_TEXTURE6);
	glBindTexture(GL_TEXTURE_2D, TextureArray[2]);

	glActiveTexture(GL_TEXTURE7);
	glBindTexture(GL_TEXTURE_2D, TextureArray[3]);

	glActiveTexture(GL_TEXTURE8);
	glBindTexture(GL_TEXTURE_2D, TextureArray[4]);

	glActiveTexture(GL_TEXTURE9);
	glBindTexture(GL_TEXTURE_2D, TextureArray[5]);

	glActiveTexture(GL_TEXTURE10);
	glBindTexture(GL_TEXTURE_2D, TextureArray[6]);

	glActiveTexture(GL_TEXTURE11);
	glBindTexture(GL_TEXTURE_CUBE_MAP, TextureArray[7]);

	glUniform1i(tex_loc, 4);
	glUniform1i(tex_loc1, 5);
	glUniform1i(tex_loc2, 6);
	glUniform1i(tex_loc3, 7);
	glUniform1i(tex_loc4, 8);
	glUniform1i(tex_loc5, 9);
	glUniform1i(tex_normalMap_loc, 10);
	glUniform1i(tex_loc_cube, 11);

	// glUniform1i(tex_sphereMap_loc, 8);
}

void aiRecursive_render(const aiScene *sc, const aiNode *nd)
{
	GLint loc;

	// Get node transformation matrix
	aiMatrix4x4 m = nd->mTransformation;
	// OpenGL matrices are column major
	m.Transpose();

	// save model matrix and apply node transformation
	pushMatrix(MODEL);

	float aux[16];
	memcpy(aux, &m, sizeof(float) * 16);
	multMatrix(MODEL, aux);

	// draw all meshes assigned to this node
	for (unsigned int n = 0; n < nd->mNumMeshes; ++n)
	{

		// send the material
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.ambient");
		glUniform4fv(loc, 1, spyMeshes[nd->mMeshes[n]].mat.ambient);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.diffuse");
		glUniform4fv(loc, 1, spyMeshes[nd->mMeshes[n]].mat.diffuse);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.specular");
		glUniform4fv(loc, 1, spyMeshes[nd->mMeshes[n]].mat.specular);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.emissive");
		glUniform4fv(loc, 1, spyMeshes[nd->mMeshes[n]].mat.emissive);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.shininess");
		glUniform1f(loc, spyMeshes[nd->mMeshes[n]].mat.shininess);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.texCount");
		glUniform1i(loc, spyMeshes[nd->mMeshes[n]].mat.texCount);

		unsigned int diffMapCount = 0; // read 2 diffuse textures

		// devido ao fragment shader suporta 2 texturas difusas simultaneas, 1 especular e 1 normal map

		glUniform1i(normalMap_loc, false); // GLSL normalMap variable initialized to 0
		glUniform1i(specularMap_loc, false);
		glUniform1ui(diffMapCount_loc, 0);

		if (spyMeshes[nd->mMeshes[n]].mat.texCount != 0)
			for (unsigned int i = 0; i < spyMeshes[nd->mMeshes[n]].mat.texCount; ++i)
			{
				if (spyMeshes[nd->mMeshes[n]].texTypes[i] == DIFFUSE)
				{
					if (diffMapCount == 0)
					{
						diffMapCount++;
						loc = glGetUniformLocation(shader.getProgramIndex(), "texUnitDiff");
						glUniform1i(loc, spyMeshes[nd->mMeshes[n]].texUnits[i]);
						glUniform1ui(diffMapCount_loc, diffMapCount);
					}
					else if (diffMapCount == 1)
					{
						diffMapCount++;
						loc = glGetUniformLocation(shader.getProgramIndex(), "texUnitDiff1");
						glUniform1i(loc, spyMeshes[nd->mMeshes[n]].texUnits[i]);

						glUniform1ui(diffMapCount_loc, diffMapCount);
					}
					else
						printf("Only supports a Material with a maximum of 2 diffuse textures\n");
				}
				else if (spyMeshes[nd->mMeshes[n]].texTypes[i] == SPECULAR)
				{
					loc = glGetUniformLocation(shader.getProgramIndex(), "texUnitSpec");
					glUniform1i(loc, spyMeshes[nd->mMeshes[n]].texUnits[i]);
					glUniform1i(specularMap_loc, true);
				}
				else if (spyMeshes[nd->mMeshes[n]].texTypes[i] == NORMALS)
				{ // Normal map
					loc = glGetUniformLocation(shader.getProgramIndex(), "texUnitNormalMap");
					if (normalMapKey)
						glUniform1i(normalMap_loc, normalMapKey);
					glUniform1i(loc, spyMeshes[nd->mMeshes[n]].texUnits[i]);
				}
				else
					printf("Texture Map not supported\n");
			}

		// send matrices to OGL
		computeDerivedMatrix(PROJ_VIEW_MODEL);
		glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
		glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
		computeNormalMatrix3x3();
		glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

		// bind VAO
		glBindVertexArray(spyMeshes[nd->mMeshes[n]].vao);

		if (!shader.isProgramValid())
		{
			printf("Program Not Valid!\n");
			exit(1);
		}
		// draw
		glDrawElements(spyMeshes[nd->mMeshes[n]].type, spyMeshes[nd->mMeshes[n]].numIndexes, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}

	// draw all children
	for (unsigned int n = 0; n < nd->mNumChildren; ++n)
	{
		aiRecursive_render(sc, nd->mChildren[n]);
	}
	popMatrix(MODEL);
}

void drawBumpedCubes()
{

	GLint loc;
	objId = 0;

	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{

			if (i == 0 || i == 3 || ((i == 1 || i == 2) && j == 2))
			{

				loc = glGetUniformLocation(shader.getProgramIndex(), "mat.ambient");
				glUniform4fv(loc, 1, recreativeMeshes[objId].mat.ambient);
				loc = glGetUniformLocation(shader.getProgramIndex(), "mat.diffuse");
				glUniform4fv(loc, 1, recreativeMeshes[objId].mat.diffuse);
				loc = glGetUniformLocation(shader.getProgramIndex(), "mat.specular");
				glUniform4fv(loc, 1, recreativeMeshes[objId].mat.specular);
				loc = glGetUniformLocation(shader.getProgramIndex(), "mat.shininess");
				glUniform1f(loc, recreativeMeshes[objId].mat.shininess);

				pushMatrix(MODEL);
				translate(MODEL, recreativeMeshes[objId].position[0], recreativeMeshes[objId].position[1] + j * recreativeMeshes[objId].size[1], recreativeMeshes[objId].position[2] + i * recreativeMeshes[objId].size[1]);
				scale(MODEL, recreativeMeshes[objId].size[0], recreativeMeshes[objId].size[1], recreativeMeshes[objId].size[2]);
				glUniform1i(texMode_uniformId, 6);

				// send matrices to OGL
				computeDerivedMatrix(PROJ_VIEW_MODEL);
				glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
				glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
				computeNormalMatrix3x3();
				glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

				// Render mesh
				glBindVertexArray(recreativeMeshes[objId].vao);
				glDrawElements(recreativeMeshes[objId].type, recreativeMeshes[objId].numIndexes, GL_UNSIGNED_INT, 0);
				glBindVertexArray(0);
				popMatrix(MODEL);
			}
		}
	}
}

void drawReflectedSpheres()
{
	pushMatrix(VIEW);

	GLint loc;
	for (int i = 1; i < recreativeMeshes.size(); i++)
	{

		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.ambient");
		glUniform4fv(loc, 1, recreativeMeshes[i].mat.ambient);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.diffuse");
		glUniform4fv(loc, 1, recreativeMeshes[i].mat.diffuse);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.specular");
		glUniform4fv(loc, 1, recreativeMeshes[i].mat.specular);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.shininess");
		glUniform1f(loc, recreativeMeshes[i].mat.shininess);
		pushMatrix(MODEL);
		translate(MODEL, recreativeMeshes[i].position[0], recreativeMeshes[i].position[1], recreativeMeshes[i].position[2]);

		glUniformMatrix4fv(view_uniformId, 1, GL_FALSE, mMatrix[VIEW]);
		computeDerivedMatrix(PROJ_VIEW_MODEL);
		glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
		glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
		computeNormalMatrix3x3();
		glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

		glUniform1i(texMode_uniformId, 7);

		glBindVertexArray(recreativeMeshes[i].vao);
		glDrawElements(recreativeMeshes[i].type, recreativeMeshes[i].numIndexes, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);

		popMatrix(MODEL);
	}
	glUniform1i(texMode_uniformId, 1);
	popMatrix(VIEW);
}

void renderFloor(bool forStencil)
{

	if (!forStencil) {
		GLint loc;
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.ambient");
		glUniform4fv(loc, 1, myMeshes[0].mat.ambient);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.diffuse");
		glUniform4fv(loc, 1, myMeshes[0].mat.diffuse);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.specular");
		glUniform4fv(loc, 1, myMeshes[0].mat.specular);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.shininess");
		glUniform1f(loc, myMeshes[0].mat.shininess);
		
		 glUniform1i(texMode_uniformId, pointlights ? 9 : 1);
	}

	pushMatrix(MODEL);
	scale(MODEL, floorX, floorY, floorZ);
	translate(MODEL, -0.5f, -0.5f, -0.5f);

	// send matrices to OGL
	computeDerivedMatrix(PROJ_VIEW_MODEL);
	glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
	glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
	computeNormalMatrix3x3();
	glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

	// Render mesh
	glBindVertexArray(myMeshes[0].vao);
	glDrawElements(myMeshes[0].type, myMeshes[0].numIndexes, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
	popMatrix(MODEL);
}

void positionEnvironmentObjects()
{

	GLint loc;

	for (int i = 1; i < myMeshes.size(); i++)
	{
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.ambient");
		glUniform4fv(loc, 1, myMeshes[i].mat.ambient);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.diffuse");
		glUniform4fv(loc, 1, myMeshes[i].mat.diffuse);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.specular");
		glUniform4fv(loc, 1, myMeshes[i].mat.specular);
		loc = glGetUniformLocation(shader.getProgramIndex(), "mat.shininess");
		glUniform1f(loc, myMeshes[i].mat.shininess);

		pushMatrix(MODEL);
		translate(MODEL, myMeshes[i].position[0], myMeshes[i].position[1], myMeshes[i].position[2]);
		scale(MODEL, myMeshes[i].size[0], myMeshes[i].size[1], myMeshes[i].size[2]);

		glUniform1i(texMode_uniformId, 1);

		// scale(MODEL, 2.0f, 3.0f, 0.5f);

		// send matrices to OGL
		computeDerivedMatrix(PROJ_VIEW_MODEL);
		glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
		glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
		computeNormalMatrix3x3();
		glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

		// Render mesh
		glBindVertexArray(myMeshes[i].vao);
		glDrawElements(myMeshes[i].type, myMeshes[i].numIndexes, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
		popMatrix(MODEL);
	}
}

void renderRearCamera()
{
	glUniform1i(texMode_uniformId, 0);
	MyMesh rearCubeCamera = createCube();

	int w = WinX;
	int h = WinY;
	/* create a diamond shaped stencil area */
	loadIdentity(PROJECTION);
	if (w <= h)
		ortho(-2.0, 2.0, -2.0 * (GLfloat)h / (GLfloat)w,
				2.0 * (GLfloat)h / (GLfloat)w, -10, 10);
	else
		ortho(-2.0 * (GLfloat)w / (GLfloat)h,
				2.0 * (GLfloat)w / (GLfloat)h, -2.0, 2.0, -10, 10);

	// load identity matrices for Model-View
	loadIdentity(VIEW);
	loadIdentity(MODEL);
	pushMatrix(MODEL);

	scale(MODEL, 3, 0.75, 1);
	translate(MODEL, -0.5f, -2.5f, -0.5f);

	// send matrices to OGL
	computeDerivedMatrix(PROJ_VIEW_MODEL);
	// glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
	glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
	computeNormalMatrix3x3();
	glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

	glBindVertexArray(rearCubeCamera.vao);
	glDrawElements(rearCubeCamera.type, rearCubeCamera.numIndexes, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

	popMatrix(MODEL);
}

void renderSkyBox(bool beforeStencil)
{
	int objId = 0;

	glUniform1i(texMode_uniformId, 3);

	glDepthMask(GL_FALSE);
	glFrontFace(GL_CW);

	pushMatrix(MODEL);
	pushMatrix(VIEW);

	mMatrix[VIEW][12] = 0.0f;
	mMatrix[VIEW][13] = 0.0f;
	mMatrix[VIEW][14] = 0.0f;

	scale(MODEL, 100.0f, 100.0f, 100.0f);
	translate(MODEL, -0.5f, -0.5f, -0.5f);

	glUniformMatrix4fv(model_uniformId, 1, GL_FALSE, mMatrix[MODEL]); // Transforma��o de modela��o do cubo unit�rio para o "Big Cube"
	computeDerivedMatrix(PROJ_VIEW_MODEL);
	glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);

	glBindVertexArray(skyBoxMesh[objId].vao);
	glDrawElements(skyBoxMesh[objId].type, skyBoxMesh[objId].numIndexes, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
	popMatrix(MODEL);
	popMatrix(VIEW);

	glFrontFace(GL_CCW); // restore counter clockwise vertex order to mean the front
	glDepthMask(GL_TRUE);
	glUniform1i(texMode_uniformId, 0);
}

void secondaryRenderScene(bool beforeStencil) {
	renderSkyBox(beforeStencil);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // Blend specular Ground with reflected geometry
	positionEnvironmentObjects();
	glDisable(GL_BLEND);

	showFlares();
}

void mainRenderScene()
{
	glUniform1i(texMode_uniformId, 2);
	pushMatrix(MODEL);
	translate(MODEL, 70.0f, 8.0f, 45.0f);
	scale(MODEL, 0.2f, 0.2f, 0.2f);
	rotate(MODEL, 180, 0, 1.0f, 0);
	aiRecursive_render(scene, scene->mRootNode);
	popMatrix(MODEL);
	glUniform1i(texMode_uniformId, 1);

	moveSceneObjects();

	drawReflectedSpheres();

	drawBumpedCubes();

	if (!shadowMode)
		drawTreeBillboards();
}

void reflection(float* dirLight, GLint dirLight_id, void (drawObjects)(), void (drawMirror)(bool))
{
	float res[4];
	// Render the reflected geometry
	dirLight[1] *= (-1.0f); // mirror the position of light
	multMatrixPoint(VIEW, dirLight, res);
	glUniform4fv(dirLight_id, 1, res);

	pushMatrix(MODEL);
	scale(MODEL, 1.0f, -1.0f, 1.0f);
	glCullFace(GL_FRONT);
	// draw_objects();
	drawObjects();
	glCullFace(GL_BACK);
	popMatrix(MODEL);

	dirLight[1] *= (-1.0f); // reset the light position
	multMatrixPoint(VIEW, dirLight, res);
	glUniform4fv(dirLight_id, 1, res);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // Blend specular Ground with reflected geometry
	drawMirror(false);
	glDisable(GL_BLEND);
}

void shadows(float* dirLight, GLfloat* plano_chao, void (*renderObjects)())
{
	shadowMode = true;
	float mat[16];
	// Render the Shadows
	glUniform1i(shadowMode_uniformId, true); // Render with constant color
	shadow_matrix(mat, plano_chao, dirLight);

	glDisable(GL_DEPTH_TEST); // To force the shadow geometry to be rendered even if behind the floor

	// Dark the color stored in color buffer
	glEnable(GL_BLEND);
	glBlendFunc(GL_DST_COLOR, GL_ZERO);
	glStencilOp(GL_KEEP, GL_KEEP, GL_ZERO);

	pushMatrix(MODEL);
	multMatrix(MODEL, mat);
	renderObjects();
	popMatrix(MODEL);
	glDisable(GL_BLEND);
	shadowMode = false;
}

void reflectionsAndShadows(GLfloat* plano_chao)
{
	reflection(pointLightPositions[3], lPos_uniformId_6, mainRenderScene, renderFloor);

	if(pointlights)
		shadows(pointLightPositions[3], plano_chao, mainRenderScene);
}

void renderScene(void)
{
	glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
	const float *roverCenter = rover.getRoverCenter();
	FrameCount++;
	glClearStencil(0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	// load identity matrices
	loadIdentity(VIEW);
	loadIdentity(MODEL);
	loadIdentity(PROJECTION);

	glEnable(GL_STENCIL_TEST);
	glEnable(GL_DEPTH_TEST);

	// use our shader
	glUseProgram(shader.getProgramIndex());
	setupTextures();
	setupLightInScene();

	if (cams[active_camera].camPos[1] > -1 && active_camera == 2) {
		setupCamera(false);

		// reflexion and shadows setup
		glStencilFunc(GL_NEVER, 0x1, 0x1);
		glStencilOp(GL_REPLACE, GL_KEEP, GL_KEEP);

		renderFloor(true);
		glUniform1i(shadowMode_uniformId, false);

		glStencilFunc(GL_EQUAL, 0x1, 0x1);
		glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

		GLfloat plano_chao[4] = { 0, 1, 0, 0 };
		reflectionsAndShadows(plano_chao);

		glDisable(GL_STENCIL_TEST);
		glDisable(GL_BLEND);
		glEnable(GL_DEPTH_TEST);

		glUniform1i(shadowMode_uniformId, false);
		mainRenderScene();
		secondaryRenderScene(false);
	}
	else if (active_camera == 3) {
		glUniform1i(shadowMode_uniformId, false);

		// rear vision
		glStencilFunc(GL_NEVER, 0x1, 0x1);
		glStencilOp(GL_REPLACE, GL_KEEP, GL_KEEP);
		renderRearCamera();
		glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
		
		glStencilFunc(GL_EQUAL, 0x1, 0x1);
		loadIdentity(VIEW); loadIdentity(PROJECTION); loadIdentity(MODEL);
		setupCamera(true);
		setupTextures();
		setupLightInScene();
		renderFloor(false);	
		mainRenderScene();
		secondaryRenderScene(true);

		glStencilFunc(GL_NOTEQUAL, 0x1, 0x1);
		loadIdentity(VIEW); loadIdentity(PROJECTION); loadIdentity(MODEL);
		setupCamera(false);
		setupTextures();
		setupLightInScene();
		renderFloor(false);
		mainRenderScene();
		secondaryRenderScene(false);
	}
	else {
		setupCamera(false);
		glDisable(GL_STENCIL_TEST);
		glUniform1i(shadowMode_uniformId, false);
		renderFloor(false);
		mainRenderScene();
		secondaryRenderScene(false);
	}


	glDisable(GL_DEPTH_TEST);
	// the glyph contains transparent background colors and non-transparent 
	// for the actual character pixels. So we use the blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	int m_viewport[4];
	glGetIntegerv(GL_VIEWPORT, m_viewport);

	pushMatrix(MODEL);
	loadIdentity(MODEL);
	pushMatrix(PROJECTION);
	loadIdentity(PROJECTION);
	pushMatrix(VIEW);
	loadIdentity(VIEW);

	stringstream ss;
	ss << points;
	string s;
	ss >> s;

	int str_length = s.length();

	for (int i = 0; i < 7 - str_length; i++)
		s = "0" + s;

	ortho(m_viewport[0], m_viewport[0] + m_viewport[2] - 1, m_viewport[1], m_viewport[1] + m_viewport[3] - 1, -1, 1);
	RenderText(shaderText, "Lives " + std::to_string(numLives), m_viewport[2] - 80, m_viewport[3] - 30, 0.5f, 0.3, 0.7f, 0.9f);
	RenderText(shaderText, s, m_viewport[2] / 180, m_viewport[3] - 30, 0.5f, 0.3, 0.7f, 0.9f);

	if (gameState == PAUSE)
	{
		// viewer at origin looking down at  negative z direction

		pushMatrix(MODEL);
		loadIdentity(MODEL);
		pushMatrix(PROJECTION);
		loadIdentity(PROJECTION);
		pushMatrix(VIEW);
		loadIdentity(VIEW);

		ortho(m_viewport[0], m_viewport[0] + m_viewport[2] - 1, m_viewport[1], m_viewport[1] + m_viewport[3] - 1, -1, 1);
		RenderText(shaderText, "Pause", m_viewport[2] / 2 - 4, m_viewport[3] / 2, 0.5f, 0.3, 0.7f, 0.9f);
	}

	else if (gameState == GAME_FINISHED)
	{
		pushMatrix(MODEL);
		loadIdentity(MODEL);
		pushMatrix(PROJECTION);
		loadIdentity(PROJECTION);
		pushMatrix(VIEW);
		loadIdentity(VIEW);

		ortho(m_viewport[0], m_viewport[0] + m_viewport[2] - 1, m_viewport[1], m_viewport[1] + m_viewport[3] - 1, -1, 1);
		RenderText(shaderText, "GAME OVER!", m_viewport[2] / 2 - 4, m_viewport[3] / 2, 0.5f, 0.3, 0.7f, 0.9f);
		RenderText(shaderText, "Press 'R' to Restart the Game", m_viewport[2] / 2 - 4, m_viewport[3] / 2 - 50, 0.5f, 0.3, 0.7f, 0.9f);
	}

	else
	{
		if (roverCenter[0] != 0 || roverCenter[2] != 0)
		{
			points += 1;
		}
	}

	if (!fireworks && points % 500 == 0 && points > 0)
	{
		iniParticles();
		fireworks = true;
	}

	popMatrix(PROJECTION);
	popMatrix(VIEW);
	popMatrix(MODEL);
	 
	// glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_STENCIL_TEST);
	// cout << scene->mRootNode<< endl;
	// glDisable(GL_BLEND);
	glutSwapBuffers();
}
// ------------------------------------------------------------
//
// Events from the Keyboard
//

void KeyboardUp(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'o':
	case 'O':
		wasUpz = true;
		movStateZ = RIGHT;
		acceleration_z = acceleration_up_z;
		break;
	case 'p':
	case 'P':
		wasUpz = true;
		movStateZ = LEFT;
		acceleration_z = acceleration_up_z;
		break;
	case 'q':
	case 'Q':
		wasUpx = true;
		movStateX = BACKWARD;
		acceleration_x = acceleration_up_x;
		break;
	case 'a':
	case 'A':
		wasUpx = true;
		movStateX = FORWARD;
		acceleration_x = acceleration_up_x;
		break;
	default:
		break;
	}
}

void processKeys(unsigned char key, int xx, int yy)
{
	if (gameState != GAME_FINISHED)
	{
		switch (key)
		{
		case 27:
			glutLeaveMainLoop();
			return;
			// case 'c':
			//	printf("Camera Spherical Coordinates (%f, %f, %f)\n", alpha, beta, r);
			//	return;
		case 'm':
			glEnable(GL_MULTISAMPLE);
			return;
		case 'k':
			glDisable(GL_MULTISAMPLE);
			return;
		case '1':
			active_camera = 0;
			billboard_type = 0;
			break;
		case '2':
			active_camera = 1;
			billboard_type = 2;
			break;
		case '3':
			active_camera = 2;
			billboard_type = 1;
			break;
		case '4':
			active_camera = 3;
			billboard_type = 1;
			break;
		case 'o':
		case 'O':
			wasUpz = false;
			movStateZ = LEFT;
			acceleration_z = acceleration_down_z;
			break;
		case 'p':
		case 'P':
			wasUpz = false;
			movStateZ = RIGHT;
			acceleration_z = acceleration_down_z;
			break;
		case 'q':
		case 'Q':
			wasUpx = false;
			movStateX = FORWARD;
			acceleration_x = acceleration_down_x;
			break;
		case 'a':
		case 'A':
			wasUpx = false;
			movStateX = BACKWARD;
			acceleration_x = acceleration_down_x;
			break;
		case 'c':
			pointlights = !pointlights;
			break;
		case 'h':
			headlights = !headlights;
			break;
		case 'n':
			sunlight = !sunlight;
			break;
		case 'f':
			fog = !fog;
			break;
		case 's':
		case 'S':
			gameState = gameState == PLAY ? TRANSITION_PLAY : TRANSITION_PAUSE;
		}
	}
	else
	{
		if (key == 'r' || key == 'R')
		{
			gameState = TRANSITION_FINISHED;
			numLives = baseNumLives;
			
		}
	}
}

// ------------------------------------------------------------
//
// Mouse Events
//

void processMouseButtons(int button, int state, int xx, int yy)
{
	// start tracking the mouse
	if (state == GLUT_DOWN)
	{
		startX = xx;
		startY = yy;
		if (button == GLUT_LEFT_BUTTON)
			tracking = 1;
		else if (button == GLUT_RIGHT_BUTTON)
			tracking = 2;
	}

	// stop tracking the mouse
	else if (state == GLUT_UP)
	{
		if (tracking == 1)
		{
			alpha -= (xx - startX);
			beta += (yy - startY);
		}
		else if (tracking == 2)
		{
			r += (yy - startY) * 0.01f;
			if (r < 0.1f)
				r = 0.1f;
		}
		tracking = 0;
	}
}

// Track mouse motion while buttons are pressed

void processMouseMotion(int xx, int yy)
{

	int deltaX, deltaY;
	float alphaAux, betaAux;
	float rAux;

	deltaX = -xx + startX;
	deltaY = yy - startY;

	// left mouse button: move camera
	if (tracking == 1)
	{

		alphaAux = alpha + deltaX;
		betaAux = beta + deltaY;

		if (betaAux > 85.0f)
			betaAux = 85.0f;
		else if (betaAux < -85.0f)
			betaAux = -85.0f;
		rAux = r;
	}
	// right mouse button: zoom
	else if (tracking == 2)
	{

		alphaAux = alpha;
		betaAux = beta;
		rAux = r + (deltaY * 0.01f);
		if (rAux < 0.1f)
			rAux = 0.1f;
	}

	if (active_camera == 2)
	{
		cams[2].camPos[0] = rAux * sin(alphaAux * 3.14f / 180.0f) * cos(betaAux * 3.14f / 180.0f);
		cams[2].camPos[2] = rAux * cos(alphaAux * 3.14f / 180.0f) * cos(betaAux * 3.14f / 180.0f);
		cams[2].camPos[1] = rAux * sin(betaAux * 3.14f / 180.0f);
	}

	//  uncomment this if not using an idle or refresh func
	//	glutPostRedisplay();
}

void mouseWheel(int wheel, int direction, int x, int y)
{

	r += direction * 0.1f;
	if (r < 0.1f)
		r = 0.1f;

	if (active_camera == 2)
	{
		cams[2].camPos[0] = r * sin(alpha * 3.14f / 180.0f) * cos(beta * 3.14f / 180.0f);
		cams[2].camPos[2] = r * cos(alpha * 3.14f / 180.0f) * cos(beta * 3.14f / 180.0f);
		cams[2].camPos[1] = r * sin(beta * 3.14f / 180.0f);
	}

	//  uncomment this if not using an idle or refresh func
	//	glutPostRedisplay();
}

// --------------------------------------------------------
//
// Shader Stuff
//

GLuint setupShaders()
{

	// Shader for models
	shader.init();
	shader.loadShader(VSShaderLib::VERTEX_SHADER, "shaders/pointlight.vert");
	shader.loadShader(VSShaderLib::FRAGMENT_SHADER, "shaders/pointlight.frag");

	// set semantics for the shader variables
	glBindFragDataLocation(shader.getProgramIndex(), 0, "colorOut");
	glBindAttribLocation(shader.getProgramIndex(), VERTEX_COORD_ATTRIB, "position");
	glBindAttribLocation(shader.getProgramIndex(), NORMAL_ATTRIB, "normal");
	glBindAttribLocation(shader.getProgramIndex(), TEXTURE_COORD_ATTRIB, "texCoord");
	glBindAttribLocation(shader.getProgramIndex(), TANGENT_ATTRIB, "tangent");
	glBindAttribLocation(shader.getProgramIndex(), BITANGENT_ATTRIB, "bitangent");

	glLinkProgram(shader.getProgramIndex());
	printf("InfoLog for Model Rendering Shader\n%s\n\n", shaderText.getAllInfoLogs().c_str());

	/*if (!shader.isProgramValid()) {
		printf("GLSL Model Program Not Valid!\n");
		exit(1);
	}*/

	texMode_uniformId = glGetUniformLocation(shader.getProgramIndex(), "texMode");
	pvm_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_pvm");
	vm_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_viewModel");
	view_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_View");
	normal_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_normal");
	model_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_Model");
	shadowMode_uniformId = glGetUniformLocation(shader.getProgramIndex(), "shadowMode");
	tex_loc = glGetUniformLocation(shader.getProgramIndex(), "texmap");
	tex_loc1 = glGetUniformLocation(shader.getProgramIndex(), "texmap1");
	tex_loc2 = glGetUniformLocation(shader.getProgramIndex(), "texmap2");
	tex_loc3 = glGetUniformLocation(shader.getProgramIndex(), "treeTexMap");
	tex_loc4 = glGetUniformLocation(shader.getProgramIndex(), "fireworksTexMap");
	tex_loc5 = glGetUniformLocation(shader.getProgramIndex(), "rollersTexMap");
	tex_loc_cube = glGetUniformLocation(shader.getProgramIndex(), "cubeMap");
	tex_loc6 = glGetUniformLocation(shader.getProgramIndex(), "texMapFlare");
	normalMap_loc = glGetUniformLocation(shader.getProgramIndex(), "normalMap");
	tex_normalMap_loc = glGetUniformLocation(shader.getProgramIndex(), "normalMapForBumping");
	/*
	tex_sphereMap_loc = glGetUniformLocation(shader.getProgramIndex(), "sphereMap");*/
	// reflect_perFragment_uniformId = glGetUniformLocation(shader.getProgramIndex(), "reflect_perFrag");

	lPos_uniformId_1 = glGetUniformLocation(shader.getProgramIndex(), "l_pos[0]");
	lPos_uniformId_2 = glGetUniformLocation(shader.getProgramIndex(), "l_pos[1]");
	lPos_uniformId_3 = glGetUniformLocation(shader.getProgramIndex(), "l_pos[2]");
	lPos_uniformId_4 = glGetUniformLocation(shader.getProgramIndex(), "l_pos[3]");
	// lPos_uniformId_5 = glGetUniformLocation(shader.getProgramIndex(), "l_pos[4]");
	// lPos_uniformId_6 = glGetUniformLocation(shader.getProgramIndex(), "l_pos[5]");
	lPos_uniformId_7 = glGetUniformLocation(shader.getProgramIndex(), "l_pos[4]");
	lPos_uniformId_8 = glGetUniformLocation(shader.getProgramIndex(), "l_pos[5]");
	lPos_uniformId_9 = glGetUniformLocation(shader.getProgramIndex(), "l_pos[6]");

	headlights_uniformId = glGetUniformLocation(shader.getProgramIndex(), "headlights_switch");
	pointlights_uniformId = glGetUniformLocation(shader.getProgramIndex(), "pointlights_switch");
	sunlight_uniformId = glGetUniformLocation(shader.getProgramIndex(), "sunlight_switch");
	fog_uniformId = glGetUniformLocation(shader.getProgramIndex(), "fog_switch");
	dir_headlight_uniformId = glGetUniformLocation(shader.getProgramIndex(), "dir_head");

	printf("InfoLog for Per Fragment Phong Lightning Shader\n%s\n\n", shader.getAllInfoLogs().c_str());

	// Shader for bitmap Text
	shaderText.init();
	shaderText.loadShader(VSShaderLib::VERTEX_SHADER, "shaders/text.vert");
	shaderText.loadShader(VSShaderLib::FRAGMENT_SHADER, "shaders/text.frag");

	glLinkProgram(shaderText.getProgramIndex());
	printf("InfoLog for Text Rendering Shader\n%s\n\n", shaderText.getAllInfoLogs().c_str());

	if (!shaderText.isProgramValid())
	{
		printf("GLSL Text Program Not Valid!\n");
		exit(1);
	}

	return (shader.isProgramLinked() && shaderText.isProgramLinked());
}

// ------------------------------------------------------------
//
// Model loading and OpenGL setup
//

int init()
{
	MyMesh amesh;
	rover = Rover();

	/* Initialization of DevIL */
	if (ilGetInteger(IL_VERSION_NUM) < IL_VERSION)
	{
		printf("wrong DevIL version \n");
		exit(0);
	}
	ilInit();

	/// Initialization of freetype library with font_name file
	freeType_init(font_name);

	// set the camera position based on its spherical coordinates
	cams[2].camPos[0] = r * sin(alpha * 3.14f / 180.0f) * cos(beta * 3.14f / 180.0f);
	cams[2].camPos[2] = r * cos(alpha * 3.14f / 180.0f) * cos(beta * 3.14f / 180.0f);
	cams[2].camPos[1] = r * sin(beta * 3.14f / 180.0f);

	glGenTextures(5, FlareTextureArray);
	Texture2D_Loader(FlareTextureArray, "crcl.tga", 0);
	Texture2D_Loader(FlareTextureArray, "flar.tga", 1);
	Texture2D_Loader(FlareTextureArray, "hxgn.tga", 2);
	Texture2D_Loader(FlareTextureArray, "ring.tga", 3);
	Texture2D_Loader(FlareTextureArray, "sun.tga", 4);

	glGenTextures(9, TextureArray);
	Texture2D_Loader(TextureArray, "sand.tif", 0);
	Texture2D_Loader(TextureArray, "checker.png", 1);
	Texture2D_Loader(TextureArray, "sand2.jpg", 2);
	Texture2D_Loader(TextureArray, "alien.png", 3);
	Texture2D_Loader(TextureArray, "particle.tga", 4);
	Texture2D_Loader(TextureArray, "stone.tga", 5);
	Texture2D_Loader(TextureArray, "normal.tga", 6);
	Texture2D_Loader(TextureArray, "spider/wal67ar_small.jpg", 22);

	const char *filenames[] = {"posx.jpg", "negx.jpg", "posy.jpg", "negy.jpg", "posz.jpg", "negz.jpg"};

	TextureCubeMap_Loader(TextureArray, filenames, 7);

	// Texture2D_Loader(TextureArray, "lightwood.tga", 2);

	float amb[] = {0.2f, 0.15f, 0.1f, 1.0f};
	float diff[] = {0.8f, 0.6f, 0.4f, 1.0f};
	float spec[] = {0.8f, 0.8f, 0.8f, 1.0f};
	float emissive[] = {0.0f, 0.0f, 0.0f, 1.0f};
	float shininess = 100.0f;
	int texcount = 0;

	float tree_spec[] = {0.2f, 0.2f, 0.2f, 1.0f};
	float tree_shininess = 10.0f;

	std::string filepath;
	strcat(model_dir, "spider");
	std::ostringstream oss;
	oss << model_dir << "/" << model_dir << ".obj";
	filepath = oss.str(); // path of OBJ file in the VS project

	strcat(model_dir, "/"); // directory path in the VS project

	// check if file exists
	ifstream fin(filepath.c_str());
	if (!fin.fail())
	{
		fin.close();
	}

	else
		printf("Couldn't open file: %s\n", filepath.c_str());

	// import 3D file into Assimp scene graph
	if (!Import3DFromFile(filepath))
		return (0);

	// creation of Mymesh array with VAO Geometry and Material
	spyMeshes = createMeshFromAssimp(scene);

	rover = Rover(&shader);
	rover.defineMatrices(
		&vm_uniformId,
		&pvm_uniformId,
		&normal_uniformId,
		&texMode_uniformId);
	rover.setY(1.25f + 2.0f * floorY);

	for (int i = 0; i < numRollers; i++)
	{
		randomRollers[i] = RandomRoller();
		// randomRollers[i].setRange(-rangeX / 2.0f, rangeX / 2.0f, -rangeZ / 2.0f, rangeZ / 2.0f);
		randomRollers[i].init(&shader, rover.getRoverCenter());
		randomRollers[i].defineMatrices(&vm_uniformId, &pvm_uniformId, &normal_uniformId, &texMode_uniformId);
		randomRollers[i].setIsPlaying(&moveRollers);
	}

	// floor
	amesh = createCube();
	float floorAmb[4] = { 0.1,0.1,0.1f,1 };
	float floorDiff[4] = { 0.8f, 0.6f, 0.4f,0.5f };
	float floorSpec[4] = { 0.1f,0.1f,0.1f,1.0f };
	memcpy(amesh.mat.ambient, floorAmb, 4 * sizeof(float));
	memcpy(amesh.mat.diffuse, floorDiff, 4 * sizeof(float));
	memcpy(amesh.mat.specular, floorSpec, 4 * sizeof(float));
	memcpy(amesh.mat.emissive, emissive, 4 * sizeof(float));
	amesh.mat.shininess = shininess;
	amesh.mat.texCount = texcount;
	amesh.position[0] = 0.0f;
	amesh.position[1] = 0.0f;
	amesh.position[2] = 0.0f;
	myMeshes.push_back(amesh);

	for (int i = 0; i < 10; i++)
	{

		amesh = createCube();
		float ambWall[] = {0.2f, 0.15f, 0.9f, 0.3f};
		float diffWall[] = {0.8f, 0.3f, 0.9f, 0.3f};
		float specWall[] = {0.8f, 0.8f, 0.9f, 0.3f};
		float emissiveWall[] = {0.0f, 0.0f, 0.0f, 0.3f};

		memcpy(amesh.mat.ambient, ambWall, 4 * sizeof(float));
		memcpy(amesh.mat.diffuse, diffWall, 4 * sizeof(float));
		memcpy(amesh.mat.specular, specWall, 4 * sizeof(float));
		memcpy(amesh.mat.emissive, emissiveWall, 4 * sizeof(float));
		amesh.mat.shininess = shininess;
		amesh.mat.texCount = texcount;

		amesh.position[0] = (rand() % (int)floorX) - floorX / 2.0f;
		amesh.position[1] = 0.0f;
		amesh.position[2] = (rand() % (int)floorZ) - floorZ / 2.0f;

		// 5.0f, 5.0f, 1.0f
		amesh.size[0] = 5.0f;
		amesh.size[1] = 5.0f;
		amesh.size[2] = 1.0f;

		myMeshes.push_back(amesh);
	}

	// amesh = createCube();
	// float ambWall[] = { 0.2f, 0.15f, 0.9f, 0.3f };
	// float diffWall[] = { 0.8f, 0.3f, 0.9f, 0.3f };
	// float specWall[] = { 0.8f, 0.8f, 0.9f, 0.3f };
	// float emissiveWall[] = { 0.0f, 0.0f, 0.0f, 0.3f };

	// memcpy(amesh.mat.ambient, ambWall, 4 * sizeof(float));
	// memcpy(amesh.mat.diffuse, diffWall, 4 * sizeof(float));
	// memcpy(amesh.mat.specular, specWall, 4 * sizeof(float));
	// memcpy(amesh.mat.emissive, emissiveWall, 4 * sizeof(float));
	// amesh.mat.shininess = shininess;
	// amesh.mat.texCount = texcount;

	// amesh.position[0] = floorX / 2.0f;
	// amesh.position[1] = 0.0f;
	// amesh.position[2] = floorZ / 2.0f;

	//// 5.0f, 5.0f, 1.0f
	// amesh.size[0] = 5.0f;
	// amesh.size[1] = 5.0f;
	// amesh.size[2] = 1.0f;

	// myMeshes.push_back(amesh);

	// create geometry and VAO of the quad for trees
	// objId= ..?;
	amesh = createQuad(6, 6);
	memcpy(amesh.mat.ambient, amb, 4 * sizeof(float));
	memcpy(amesh.mat.diffuse, diff, 4 * sizeof(float));
	memcpy(amesh.mat.specular, tree_spec, 4 * sizeof(float));
	memcpy(amesh.mat.emissive, emissive, 4 * sizeof(float));
	amesh.mat.shininess = tree_shininess;
	amesh.mat.texCount = texcount;
	billboardMesh.push_back(amesh);

	// create geometry and VAO of the quad for particles
	// objId = ..?;
	amesh = createQuad(2, 2);
	amesh.mat.texCount = texcount;
	billboardMesh.push_back(amesh);

	amesh = createCube();
	memcpy(amesh.mat.ambient, amb, 4 * sizeof(float));
	memcpy(amesh.mat.diffuse, diff, 4 * sizeof(float));
	memcpy(amesh.mat.specular, spec, 4 * sizeof(float));
	memcpy(amesh.mat.emissive, emissive, 4 * sizeof(float));
	amesh.mat.shininess = shininess;
	amesh.mat.texCount = texcount;
	skyBoxMesh.push_back(amesh);

	rearCubeCamera = createCube();
	memcpy(rearCubeCamera.mat.ambient, amb, 4 * sizeof(float));
	memcpy(rearCubeCamera.mat.diffuse, diff, 4 * sizeof(float));
	memcpy(rearCubeCamera.mat.specular, spec, 4 * sizeof(float));
	memcpy(rearCubeCamera.mat.emissive, emissive, 4 * sizeof(float));
	rearCubeCamera.mat.shininess = shininess;
	rearCubeCamera.mat.texCount = texcount;

	amesh = createCube();
	memcpy(amesh.mat.ambient, amb, 4 * sizeof(float));
	memcpy(amesh.mat.diffuse, diff, 4 * sizeof(float));
	memcpy(amesh.mat.specular, spec, 4 * sizeof(float));
	memcpy(amesh.mat.emissive, emissive, 4 * sizeof(float));
	amesh.mat.shininess = shininess;
	amesh.mat.texCount = texcount;
	amesh.position[0] = 8.0f;
	amesh.position[1] = 0.0f;
	amesh.position[2] = 10.0f;
	amesh.size[0] = 3.0f;
	amesh.size[1] = 5.0f;
	amesh.size[2] = 5.0f;
	recreativeMeshes.push_back(amesh);

	for (int i = 0; i < 10; i++)
	{

		amesh = createSphere(2.0f, 20);
		memcpy(amesh.mat.ambient, amb, 4 * sizeof(float));
		memcpy(amesh.mat.diffuse, diff, 4 * sizeof(float));
		memcpy(amesh.mat.specular, spec, 4 * sizeof(float));
		memcpy(amesh.mat.emissive, emissive, 4 * sizeof(float));
		amesh.mat.shininess = shininess;
		amesh.mat.texCount = texcount;
		amesh.position[0] = (rand() % (int)floorX) - floorX / 2.0f;
		amesh.position[1] = 2.0f;
		amesh.position[2] = (rand() % (int)floorZ) - floorZ / 2.0f;
		recreativeMeshes.push_back(amesh);
	}

	// create geometry and VAO of the quad for flare elements
	amesh = createQuad(1, 1);
	flareMeshes.push_back(amesh);

	// Load flare from file
	loadFlareFile(&AVTflare, "flare.txt");

	// some GL settings
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearStencil(0x0);
	glEnable(GL_STENCIL_TEST);
	glEnable(GL_MULTISAMPLE);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	// Codigo do prof faz:
	/* glEnable(GL_CULL_FACE);
	glEnable(GL_MULTISAMPLE);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f); */

	// top perspective pos: (0,20,0)
	cams[0].camPos[1] = 20.0f;

	// top ortho pos: (0,20,0)
	cams[1].camPos[1] = 20.0f;
	cams[1].mode = 1;
	return (1);
}

int sign(float x) { return (x > 0.0f) ? 1 : -1; }

void fixedRenderScene(int value)
{

	rover.move(velocity_x, velocity_z, 1.0f / ((float) value));
	switch (gameState)
	{
	case TRANSITION_PLAY:
		pause_velocity_x = velocity_x;
		pause_velocity_z = velocity_z;
		velocity_x = 0;
		velocity_z = 0;
		gameState = PAUSE;
		moveRollers = false;
		break;

	case TRANSITION_PAUSE:
		velocity_x = pause_velocity_x;
		velocity_z = pause_velocity_z;
		gameState = PLAY;
		moveRollers = true;
		break;

	case TRANSITION_FINISHED:
		points = 0;
		velocity_x = 0;
		velocity_z = 0;
		acceleration_z = acceleration_down_z;
		acceleration_x = acceleration_down_x;
		movStateX = IDLE;
		movStateZ = IDLE;
		gameState = PLAY;
		wasUpx = false;
		wasUpz = false;
		break;

	case PLAY:
		// se as teclas q ou a n�o est�o a ser carregadas,
		// leva a vel_x at� zero

		if (wasUpx)
		{
			if (abs(velocity_x) < velocity_x)
			{
				velocity_x = 0.0f;
				movStateX = IDLE;
				wasUpx = false;
			}
			else
			{
				velocity_x += -sign(velocity_x) * acceleration_x * (1.0f / ((float)value));
			}
		}
		else if (abs(velocity_x) >= velocity_max_x)
		{
			velocity_x = velocity_x > 0 ? velocity_max_x * 0.95f : -velocity_max_x * 0.95f;
		}
		else
		{
			switch (movStateX)
			{
			case FORWARD:
				velocity_x += acceleration_x * (1.0f / ((float)value));
				break;
			case BACKWARD:
				velocity_x -= acceleration_x * (1.0f / ((float)value));
				break;
			default:
				break;
			}
		}

		if (wasUpz)
		{
			if (abs(velocity_z) < velocity_min)
			{
				velocity_z = 0.0f;
				movStateZ = IDLE;
				wasUpz = false;
			}
			else if (abs(velocity_z) >= velocity_max_z)
			{
				velocity_z = velocity_z > 0 ? velocity_max_z * 0.95f : -velocity_max_z * 0.95f;
			}
			else
			{
				velocity_z += -sign(velocity_z) * acceleration_z * (1.0f / ((float)value));
			}
		}
		else
		{
			switch (movStateZ)
			{
			case LEFT:
				velocity_z += acceleration_z * (1.0f / ((float)value));
				break;
			case RIGHT:
				velocity_z -= acceleration_z * (1.0f / ((float)value));
				break;
			default:
				break;
			}
		}

		for (int i = 0; i < numRollers; i++)
		{
			randomRollers[i].incrementVelocity(1.0f / value);
		}
		break;

	case GAME_FINISHED:
	{
		rover.resetRover();
	}
	default:
		break;
	}

	const float *roverPosition = rover.getRoverCenter();
	test_file << roverPosition[0] << "," << roverPosition[1] << "," << roverPosition[2] << "," << velocity_x << "," << velocity_z << "," << deltaT << "," << endl;

	for (int i = 0; i < numRollers; i++)
	{
		randomRollers[i].moveRoller(deltaT, rover.getRoverCenter());
	}

	glutTimerFunc(1000 / value, fixedRenderScene, value);
}

void collisionDetection(int value)
{
	const float *roverCenter = rover.getRoverCenter();
	// for (int i = 0; i < numRollers; i++) {
	//	cout << "Position x " << randomRollers[i].pos[0] << endl;
	//	cout << "Position y " << randomRollers[i].pos[1] << endl;
	//	cout << "Position z " << randomRollers[i].pos[2] << endl;
	//	cout << "Radius" << randomRollers[i].radius << endl;
	// }
	float posX = roverCenter[0] - (rover.bodyX / 2);
	float posZ = roverCenter[2] - (rover.bodyZ / 2);

	for (int i = 0; i < numRollers; i++)
	{

		bool collisionX = posX + rover.bodyX >= randomRollers[i].pos[0] &&
						  randomRollers[i].pos[0] + randomRollers[i].radius * 2 >= posX;

		bool collisionZ = posZ + rover.bodyZ >= randomRollers[i].pos[2] &&
						  randomRollers[i].pos[2] + randomRollers[i].radius * 2 >= posZ;

		if (collisionX && collisionZ)
		{
			cout << "Collision Detected: Rover & Roller" << endl;
			if (roverCenter[0] != 0 || roverCenter[2] != 0)
			{
				numLives -= 1;

				if (numLives == 0)
				{

					gameState = GAME_FINISHED;
					fireworks = false;
				}
			}
			else
			{
				rover.resetRover();
				velocity_x = 0;
				velocity_z = 0;
				acceleration_z = 10.0f;
				acceleration_x = 1.0f;
				movStateX = IDLE;
				movStateZ = IDLE;
			}
		}
	}

	for (int i = 1; i < myMeshes.size(); i++)
	{

		bool collisionX = posX + rover.bodyX >= myMeshes[i].position[0] &&
						  myMeshes[i].position[0] + myMeshes[i].size[0] >= posX;

		bool collisionZ = posZ + rover.bodyZ >= myMeshes[i].position[2] &&
						  myMeshes[i].position[2] + 2 * myMeshes[i].size[2] >= posZ;

		if (collisionX && collisionZ)
		{
			cout << "hererer" << endl;
			myMeshes[i].position[0] += velocity_x * cosf(*rover.getRoverAngle() * 3.14159265358979323846f / 180.0f);
			myMeshes[i].position[2] += velocity_x * -sinf(*rover.getRoverAngle() * 3.14159265358979323846f / 180.0f);

			velocity_x = velocity_x * 0;
			velocity_z = velocity_z * 0;
			acceleration_z = 10.0f;
			acceleration_x = 1.0f;
			movStateX = IDLE;
			movStateZ = IDLE;
		}

		// if (collisionX && collisionZ) {
		//	cout << "Collided with wall" << endl;
		//	cout << "Positions " << myMeshes[i].position[0] << " " << posX << " " << myMeshes[i].position[2] << " " << posZ << endl;

		//	if (movStateX == FORWARD)
		//	{
		//		movStateX = BACKWARD;
		//	}

		//	else {
		//		movStateX = FORWARD;
		//	}

		//	if (movStateZ == LEFT)
		//	{
		//		movStateX = RIGHT;
		//	}

		//	else {
		//		movStateX = LEFT;

		//	}

		//	//velocity_x = -velocity_x * .1f;
		//	//velocity_z = -velocity_z * .1f;

		//	rover.move(velocity_x, velocity_z, deltaT);

		//	myMeshes[i].position[0] += velocity_x * cosf(*rover.getRoverAngle() * 3.14159265358979323846f / 180.0f) * deltaT;
		//	myMeshes[i].position[2] += velocity_x * -sinf(*rover.getRoverAngle() * 3.14159265358979323846f / 180.0f) * deltaT;
	}

	glutTimerFunc(1, collisionDetection, 0);
}

// ------------------------------------------------------------
//
// Main function
//

int main(int argc, char **argv)
{
	srand(time(0));

	//  GLUT initialization
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_MULTISAMPLE | GLUT_STENCIL);

	glutInitContextVersion(4, 3);
	glutInitContextProfile(GLUT_CORE_PROFILE);
	glutInitContextFlags(GLUT_FORWARD_COMPATIBLE | GLUT_DEBUG);

	glutInitWindowPosition(100, 10);
	glutInitWindowSize(WinX, WinY);
	WindowHandle = glutCreateWindow(CAPTION);

	//  Callback Registration
	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);

	glutTimerFunc(0, timer, 0);
	// glutIdleFunc(renderScene);  // Use it for maximum FPS (no cap)
	glutTimerFunc(0, refresh, 120); // use it to to get 60 FPS whatever
	glutTimerFunc(0, fixedRenderScene, 100);
	glutTimerFunc(0, collisionDetection, 0);

	//	Mouse and Keyboard Callbacks
	// glutSetKeyRepeat(GLUT_KEY_REPEAT_OFF);
	glutKeyboardFunc(processKeys);
	glutKeyboardUpFunc(KeyboardUp);
	glutMouseFunc(processMouseButtons);
	glutMotionFunc(processMouseMotion);
	glutMouseWheelFunc(mouseWheel);

	//	return from main loop
	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);

	//	Init GLEW
	glewExperimental = GL_TRUE;
	glewInit();

	printf("Vendor: %s\n", glGetString(GL_VENDOR));
	printf("Renderer: %s\n", glGetString(GL_RENDERER));
	printf("Version: %s\n", glGetString(GL_VERSION));
	printf("GLSL: %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));

	if (!setupShaders())
		return (1);

	init();

	//  GLUT main loop
	glutMainLoop();

	test_file.close();

	return (0);
}
