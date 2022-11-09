#include "AVTmathLib.h"
#include "vsShaderLib.h"
#include "geometry.h"
#include "Rover.h"



#ifdef _WIN32
#define M_PI       3.14159265358979323846f
#endif

// The storage for matrices
extern float mCompMatrix[COUNT_COMPUTED_MATRICES][16];

// The normal matrix
extern float mNormal3x3[9];

Rover::Rover() {}

Rover::Rover(VSShaderLib* shader) {
	this->shader = shader;
	this->rover[0] = createSphere(wRadius/10.0f, wSides);
	setupMesh(&rover[0]);
	this->rover[1] = createCube();
	setupMesh(&rover[1]);
	this->rover[2] = createCylinder(wHeight, wRadius, wSides);
	setupMesh(&rover[2]);
	this->rover[3] = createCylinder(wHeight, wRadius, wSides);
	setupMesh(&rover[3]);
	this->rover[4] = createCylinder(wHeight, wRadius, wSides);
	setupMesh(&rover[4]);
	this->rover[5] = createCylinder(wHeight, wRadius, wSides);
	setupMesh(&rover[5]);
	this->rover[6] = createCube();
	setupMesh(&rover[6]);
}

void Rover::setupMesh(MyMesh* mesh) {
	float amb[] = { 0.2f, 0.15f, 0.1f, 1.0f };
	float diff[] = { 0.8f, 0.6f, 0.4f, 1.0f };
	float spec[] = { 0.8f, 0.8f, 0.8f, 1.0f };
	float emissive[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	float shininess = 100.0f;
	int texcount = 0;
	memcpy(mesh->mat.ambient, amb, 4 * sizeof(float));
	memcpy(mesh->mat.diffuse, diff, 4 * sizeof(float));
	memcpy(mesh->mat.specular, spec, 4 * sizeof(float));
	memcpy(mesh->mat.emissive, emissive, 4 * sizeof(float));
	mesh->mat.shininess = shininess;
	mesh->mat.texCount = texcount;
}

void Rover::renderMesh(MyMesh* mesh) {
	GLint loc;
	loc = glGetUniformLocation(shader->getProgramIndex(), "mat.ambient");
	glUniform4fv(loc, 1, mesh->mat.ambient);
	loc = glGetUniformLocation(shader->getProgramIndex(), "mat.diffuse");
	glUniform4fv(loc, 1, mesh->mat.diffuse);
	loc = glGetUniformLocation(shader->getProgramIndex(), "mat.specular");
	glUniform4fv(loc, 1, mesh->mat.specular);
	loc = glGetUniformLocation(shader->getProgramIndex(), "mat.shininess");
	glUniform1f(loc, mesh->mat.shininess);

	
}

void Rover::defineMatrices(GLint* vmId, GLint* pvmId, GLint* normalId,GLint* texMode_uniformId) {
	this->vmId = vmId;
	this->pvmId = pvmId;
	this->normalId = normalId;
	this->tex = texMode_uniformId;
}

void Rover::sendMatrices() const {
	computeDerivedMatrix(PROJ_VIEW_MODEL);
	glUniformMatrix4fv(*vmId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
	glUniformMatrix4fv(*pvmId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
	computeNormalMatrix3x3();
	glUniformMatrix3fv(*normalId, 1, GL_FALSE, mNormal3x3);
	glUniform1i(*tex, 1);

}

void Rover::setMinimumVelocity(float minVel) { this->minVel = minVel; }

void Rover::setY(float new_y) { this->roverCenter[1] = new_y; }

//const float* Rover::getHeadLighter() const {
//	return headLighter;
//}
//
//const float* Rover::getHeadLighterV() const {
//	return headLighterV;
//}


void Rover::renderRover() {
	pushMatrix(MODEL);

	translate(MODEL, roverCenter[0], roverCenter[1], roverCenter[2]);
	rotate(MODEL, roverAngle, 0.0f, 1.0f, 0.0f);

	// build rover
	MyMesh m = rover[0]; // the rover's CG 
	renderMesh(&m);

	sendMatrices();

	glBindVertexArray(m.vao);

	glDrawElements(m.type, m.numIndexes, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

	for (int i = 1; i < 7; i++)
	{
		m = rover[i];
		renderMesh(&m);
		pushMatrix(MODEL);

		switch (i)
		{
		case 1:
			translate(MODEL, -bodyX / 2.0f, -bodyY / 2.0f, -bodyZ / 2.0f);
			scale(MODEL, bodyX, bodyY, bodyZ);
			break;
		case 2:
			translate(
				MODEL, 
				bodyX / 2.0f - wRadius - 0.1f, 
				-bodyY / 2.0f,
				bodyZ / 2.0f + wHeight / 2.0f
			);
			rotate(MODEL, 90.0f, 1.0f, 0.0f, 0.0f);
			break;
		case 3:
			translate(
				MODEL,
				bodyX / 2.0f - wRadius - 0.1f,
				-bodyY / 2.0f,
				-bodyZ / 2.0f - wHeight / 2.0f
			);
			rotate(MODEL, 90.0f, 1.0f, 0.0f, 0.0f);
			break;
		case 4:
			translate(
				MODEL,
				-bodyX / 2.0f + wRadius + 0.1f,
				-bodyY / 2.0f,
				bodyZ / 2.0f + wHeight / 2.0f
			);
			rotate(MODEL, 90.0f, 1.0f, 0.0f, 0.0f);
			break;
		case 5:
			translate(
				MODEL,
				-bodyX / 2.0f + wRadius + 0.1f,
				-bodyY / 2.0f,
				-bodyZ / 2.0f - wHeight / 2.0f
			);
			rotate(MODEL, 90.0f, 1.0f, 0.0f, 0.0f);
			break;
		case 6:
			translate(
				MODEL, 
				bodyX / 2.0f - headX - 0.1f, 
				(bodyY - headY) / 2.0f, 
				-headZ / 2.0f
			);
			scale(MODEL, headX, headY, headZ);
			break;
		default:
			break;
		}

		sendMatrices();
		glBindVertexArray(m.vao);
		glDrawElements(m.type, m.numIndexes, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
		popMatrix(MODEL);
	}

	popMatrix(MODEL);
}

void Rover::move(float velX, float velZ, float dt) {
	roverCenter[0] += velX * cosf(roverAngle * M_PI / 180.0f) * dt;
	roverCenter[2] += velX * -sinf(roverAngle * M_PI / 180.0f) * dt;
	roverAngle += velZ * dt;
}

void Rover::setInitialPosition(float posx, float posy, float posz) {
	this->roverCenter[0] = posx;
	this->roverCenter[1] = posy;
	this->roverCenter[2] = posz;
}

const float* Rover::getRoverCenter() const {
	return roverCenter;
}

const float* Rover::getRoverAngle() const {
	return &roverAngle;
}

void Rover::resetRover() {
	translate(MODEL, -roverCenter[0], -roverCenter[1], -roverCenter[2]);
	roverCenter[0] = 0;
	roverCenter[2] = 0;
	rotate(MODEL, -roverAngle, 0.0f, 1.0f, 0.0f);
	roverAngle = 0;
	

	// build rover
}

static void multMatrixPoint(float* matrix, float* point, float* res) {

	for (int i = 0; i < 4; ++i) {

		res[i] = 0.0f;

		for (int j = 0; j < 4; j++) {

			res[i] += point[j] * matrix[j + i * 4];
		}
	}
}

const float* Rover::getRearCameraConfiguration() {
	float radRoverAngle = roverAngle * M_PI / 180.0f;
	float rearCameraVectorPos[4] = {
		- bodyX / 2.0f - 0.5f,
		- bodyY / 2.0f,
		0,
		1
	};
	float rearCameraVectorAt[4] = { -10, 4.5, 0, 0 };
	float rover2WMatrix[16] = {
		cosf(radRoverAngle), 0, sinf(radRoverAngle), roverCenter[0],
						  0, 1, 				  0, roverCenter[1],
	   -sinf(radRoverAngle), 0, cosf(radRoverAngle), roverCenter[2],
						  0, 0,					  0,			  1
	};

	// position of the rear camera
	float resPos[4];
	multMatrixPoint(rover2WMatrix, rearCameraVectorPos, resPos);
	
	rearCameraConfig[0]  = resPos[0];
	rearCameraConfig[1]  = resPos[1];
	rearCameraConfig[2]  = resPos[2];

	// at of the rear camera
	float resVec[4];
	multMatrixPoint(rover2WMatrix, rearCameraVectorAt, resVec);
	rearCameraConfig[3]  = resPos[0] + resVec[0];
	rearCameraConfig[4]  = resPos[1] + resVec[1];
	rearCameraConfig[5]  = resPos[2] + resVec[2];

	// position of the front camera
	rearCameraVectorPos[0] *= -1.0f;
	multMatrixPoint(rover2WMatrix, rearCameraVectorPos, resPos);
	rearCameraConfig[6]  = resPos[0];
	rearCameraConfig[7]  = resPos[1];
	rearCameraConfig[8]  = resPos[2];

	// at of the front camera	
	rearCameraVectorAt[0] *= -1.0f;
	rearCameraVectorAt[1] = 0;
	multMatrixPoint(rover2WMatrix, rearCameraVectorAt, resVec);
	rearCameraConfig[9]  = resPos[0] + resVec[0];
	rearCameraConfig[10] = resPos[1] + resVec[1];
	rearCameraConfig[11] = resPos[2] + resVec[2];

	return rearCameraConfig;
}

const float* Rover::getHeadLightConfiguration() {
	float radRoverAngle = (roverAngle - 90) * M_PI / 180.0f;
	float rover2WMatrix[16] = {
		cosf(radRoverAngle), 0, sinf(radRoverAngle), roverCenter[0],
						  0, 1, 				  0, roverCenter[1],
	   -sinf(radRoverAngle), 0, cosf(radRoverAngle), roverCenter[2],
						  0, 0,					  0,			  1
	};
	float headLightPos[4] = {
		bodyX / 2.0f,
		bodyY / 2.0f,
		bodyZ / 2.0f,
		1.0f
	};
	float headLightAt[4] = {
		10,
		0,
		0,
		0
	};

	// left headlight
	float res[4];
	multMatrixPoint(rover2WMatrix, headLightPos, res);
	headLightConfig[0] = res[0];
	headLightConfig[1] = res[1];
	headLightConfig[2] = res[2];

	// right headLight
	headLightConfig[2] *= -1.0f;
	multMatrixPoint(rover2WMatrix, headLightPos, res);
	headLightConfig[3] = res[0];
	headLightConfig[4] = res[1];
	headLightConfig[5] = res[2];

	// direction of the headlights
	multMatrixPoint(rover2WMatrix, headLightAt, res);
	headLightConfig[6] = res[0];
	headLightConfig[7] = res[1];
	headLightConfig[8] = res[2];

	return headLightConfig;
}




