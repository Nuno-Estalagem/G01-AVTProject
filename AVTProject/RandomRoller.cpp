#include <stdlib.h>
#include <time.h>
#include "AVTmathLib.h"
#include "vsShaderLib.h"
#include "geometry.h"
#include "RandomRoller.h"

// The storage for matrices
extern float mCompMatrix[COUNT_COMPUTED_MATRICES][16];

// The normal matrix
extern float mNormal3x3[9];

#ifdef _WIN32
#define M_PI       3.14159265358979323846f
#endif

RandomRoller::RandomRoller() {}

void RandomRoller::setY(float posY) {
	this->pos[1] = posY;
}

void RandomRoller::setRadiusRange(float minR, float maxR) {
	this->limitsR[0] = minR;
	this->limitsR[1] = maxR;
}

void RandomRoller::setVelocityInternal(int mode) {
	float maxV = 0.1f;
	switch (mode) {
	case 0: vel[0] = -maxV + 2.0f * randf() * maxV; break;
	case 2: vel[2] = -maxV + 2.0f * randf() * maxV; break;
	case 4: setVelocityInternal(0); setVelocityInternal(2); break;
	default: break;
	}
}

void RandomRoller::defineMeshRoller() {
	float amb[] = { 0.2f, 0.15f, 0.1f, 1.0f };
	float diff[] = { 0.8f, 0.6f, 0.4f, 1.0f };
	float spec[] = { 0.8f, 0.8f, 0.8f, 1.0f };
	float emissive[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	float shininess = 500.0f;
	int texcount = 0;
	memcpy(this->roller.mat.ambient, amb, 4 * sizeof(float));
	memcpy(this->roller.mat.diffuse, diff, 4 * sizeof(float));
	memcpy(this->roller.mat.specular, spec, 4 * sizeof(float));
	memcpy(this->roller.mat.emissive, emissive, 4 * sizeof(float));
	this->roller.mat.shininess = shininess;
	this->roller.mat.texCount = texcount;
}

void RandomRoller::init(VSShaderLib* shader, const float* roverCenter) {
	this->shader = shader;
	this->radius = this->limitsR[0] + randf() * (this->limitsR[1] - this->limitsR[0]);
	this->roller = createSphere(this->radius, this->divisions);
	defineMeshRoller(); 
	
	// set position
	float theta = 2.0f * M_PI * randf();
	float a[3] = {
		roverCenter[0] + roverRadius * cosf(theta),
		radius,
		roverCenter[2] + roverRadius * sinf(theta)
	};
	memcpy(pos, a, 3 * sizeof(float));

	//set velocity b(roverCenter -> randomPointOnCircle) = a(roverCenter -> roller) + v(velocidade) <=> v = b - a
	float b[3] = {
		(0.1f + randf() * 0.9f) * roverRadius,
		0.0f,
		(0.1f + randf() * 0.9f) * roverRadius
	};
	normalize(b);
	float v[3] = { b[0] - a[0], 0.0f, b[2] - a[2] };
	normalize(v);
	vel[0] = globalVel * v[0];
	vel[1] = 0.0f;
	vel[2] = globalVel * v[2];
}

void RandomRoller::renderRoller() {
	GLint loc;
	loc = glGetUniformLocation(this->shader->getProgramIndex(), "mat.ambient");
	glUniform4fv(loc, 1, this->roller.mat.ambient);
	loc = glGetUniformLocation(shader->getProgramIndex(), "mat.diffuse");
	glUniform4fv(loc, 1, this->roller.mat.diffuse);
	loc = glGetUniformLocation(shader->getProgramIndex(), "mat.specular");
	glUniform4fv(loc, 1, this->roller.mat.specular);
	loc = glGetUniformLocation(shader->getProgramIndex(), "mat.shininess");
	glUniform1f(loc, this->roller.mat.shininess);


}

void RandomRoller::defineMatrices(GLint* vmId, GLint* pvmId, GLint* normalId, GLint* texMode_uniformId) {
	this->vmId = vmId;
	this->pvmId = pvmId;
	this->normalId = normalId;
	this->tex = texMode_uniformId;

}

void RandomRoller::sendMatrices() const {
	computeDerivedMatrix(PROJ_VIEW_MODEL);
	glUniformMatrix4fv(*vmId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
	glUniformMatrix4fv(*pvmId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
	computeNormalMatrix3x3();
	glUniformMatrix3fv(*normalId, 1, GL_FALSE, mNormal3x3);
	glUniform1i(*tex, 6);

}

void RandomRoller::setIsPlaying(bool* isPlaying) {
	this->isPlaying = isPlaying;
}

void RandomRoller::setPositionInternal(int axis) {
	switch (axis) {
	case 0: this->pos[0] = limitsX[0] + randf() * (limitsX[1] - limitsX[0]); break;
	case 2: this->pos[2] = limitsZ[0] + randf() * (limitsZ[1] - limitsZ[0]); break;
	case 4: 
		this->setPositionInternal(0); 
		this->setPositionInternal(2); 
		break;
	default: break;
	}
}

void RandomRoller::setRange(float minX, float maxX, float minZ, float maxZ) {
	this->limitsX[0] = minX;
	this->limitsX[1] = maxX;
	this->limitsZ[0] = minZ;
	this->limitsZ[1] = maxZ;
	this->setPositionInternal(4);
}

float RandomRoller::randf() {
	return ((float)rand()) / ((float)RAND_MAX);
}

float absv(float* a) {
	return sqrtf(
		powf(a[0], 2.0f) + 
		powf(a[1], 2.0f) + 
		powf(a[2], 2.0f)
	);
}

float distv(float* a, float* b) {
	return sqrtf(
		(a[0] - b[0]) * (a[0] - b[0]) +
		(a[1] - b[1]) * (a[1] - b[1]) +
		(a[2] - b[2]) * (a[2] - b[2])
	);
}

void RandomRoller::moveRoller(float dt, const float* roverCenter) {
	float currentDistance = distv((float*) roverCenter, pos);

	// outside of the rover circle
	if (currentDistance > roverRadius) {
		float theta = 2.0f * M_PI * randf();
		float a[3] = {
			roverCenter[0] + roverRadius * cosf(theta),
			radius,
			roverCenter[2] + roverRadius * sinf(theta)
		};
		memcpy(pos, a, 3 * sizeof(float));
		float b[3] = {
			(0.1f + randf() * 0.9f) * roverRadius,
			0.0f,
			(0.1f + randf() * 0.9f) * roverRadius
		};
		float v[3] = { b[0] - a[0], 0.0f, b[2] - a[2] };
		normalize(v);
		vel[0] = globalVel * v[0];
		vel[1] = 0.0f;
		vel[2] = globalVel * v[2];
	}

	bool isMoving = *isPlaying;
	this->pos[0] += vel[0] * dt * isMoving;
	this->pos[2] += vel[2] * dt * isMoving;

	crossProduct(vel, up, axis); normalize(axis);

	angle -= isMoving * 6.0f * absv(vel) * dt / radius;
}

void RandomRoller::drawRoller() {
	renderRoller();

	pushMatrix(MODEL);

	translate(MODEL, this->pos[0], this->pos[1], this->pos[2]);
	rotate(MODEL, angle, axis[0], axis[1], axis[2]);

	sendMatrices();

	glBindVertexArray(this->roller.vao);

	glDrawElements(this->roller.type, this->roller.numIndexes, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

	popMatrix(MODEL);
}

void RandomRoller::incrementVelocity(float dt) {
	if (increment > incdT) {
		globalVel *= (1.0f + 0.1f * (absv(vel) < maximumVel));
		increment = 0.0f;
	} else {
		increment += dt;
	}
}

const float* RandomRoller::getPosition() const {
	return this->pos;
}
