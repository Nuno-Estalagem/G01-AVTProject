#include "geometry.h"
#include <GL/glew.h>

class MyMesh;
class VSShaderLib;

#pragma once
class RandomRoller {
private:
	MyMesh roller;
	VSShaderLib* shader;
	float roverRadius = 50.0f;
	float limitsR[2] = { 0.5f, 2.5f };
	float limitsX[2] = { 0.0f, 0.0f };
	float limitsZ[2] = { 0.0f, 0.0f };
	float up[3] = { 0.0f,1.0f,0.0f };
	float angle = 0.0f;
	float axis[3] = { 0.0f };

	int divisions = 20;
	GLint* vmId;
	GLint* pvmId;
	GLint* normalId;
	GLint* tex;
	bool doesRoll = false;
	float vel[3] = { 0.0f, 0.0f, 0.0f };
	float globalVel = 10.0f;
	float maximumVel = 10.0f;
	float incdT = 1.0f;
	float increment = 0;

	bool* isPlaying;

	// Outputs a random float between 0.0f and 1.0f
	float randf();

	void renderRoller();
	void sendMatrices() const;
	void setPositionInternal(int axis);
	void setVelocityInternal(int axis);
	void defineMeshRoller();
public:
	RandomRoller();
	float radius = 0.0f;
	float pos[3] = { 0.0f, 1.5f, 0.0f };
	void incrementVelocity(float dt);
	void init(VSShaderLib* shader, const float* roverCenter);
	void setRange(float minX, float maxX, float minZ, float maxZ);
	const float* getPosition() const;
	void setRadiusRange(float minR, float maxR);
	void setY(float posY);
	void defineMatrices(GLint* vmId, GLint* pvmId, GLint* normalId, GLint* uniform_id);
	void drawRoller();
	void setIsPlaying(bool* isPlaying);
	void moveRoller(float dt, const float* actualRoverCenter);
};