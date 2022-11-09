
class MyMesh;
class VSShaderLib;

#pragma once
class Rover
{
public:
	Rover();
	Rover(VSShaderLib* shader);
	void renderMesh(MyMesh* mesh);
	void drawMesh(MyMesh* mesh);

	// The move function of the rover, where given a velocity vector
	// in terms of (x,y,z), an angle velocity in the direction of the 
	// y-axis and the time interval (deltaT) that those velocities
	// were applied, it makes the rover move through the screen. 
	// 
	// The latter may be 1000 / FPS seconds.
	void move(float velX, float velAngle, float dt);
	void renderRover();
	
	void defineMatrices(GLint* vmId, GLint* pvmId, GLint* normalId, GLint* uniform_id);

	// Get the rover's current center vector in style of (x,y,z)
	const float* getRoverCenter() const;

	// Get the rover's current angle value in degrees
	const float* getRoverAngle() const;
	//const float* getHeadLighter() const;
	//const float* getHeadLighterV() const;
	void setY(float new_y);
	void setInitialPosition(float posx, float posy, float posz);
	void setMinimumVelocity(float minVel);
	void resetRover();
	const float* getRearCameraConfiguration();
	const float* getHeadLightConfiguration();


	// The body scale for the x-axis 
	float bodyX = 4.5f;

	// The body scale for the y-axis 
	float bodyY = 1.0f;

	// The body scale for the z-axis 
	float bodyZ = 2.0f;

private:
	MyMesh rover[7];
	VSShaderLib* shader;
	GLint* vmId;
	GLint* pvmId;
	GLint* normalId;
	GLint* tex;
	
	// The wheels radius
	float wRadius = 0.75f;
	
	// The wheels height
	float wHeight = 0.25f;
	
	// The wheels side
	int wSides = 20;
	
	// The head scale for the x-axis 
	float headX = 0.5f;
	
	// The head scale for the y-axis 
	float headY = 1.5f;

	// The head scale for the z-axis 
	float headZ = 2.0f;

	// The center of the rover, that is, the center 
	// of the body
	float roverCenter[3] = { 0.0f, 0.0f, 0.0f };

	// The headlights' centers
	float headLighter[6] = { 0.0f };

	// The headlights' directions
	float headLighterV[4] = { 0.0f };

	float headLightConfig[9] = { 0 };

	// Height of the rover
	float roverHeight = 1.55f;

	// The angle between the x-axis of the world coordinates
	// and the rover x-axis. This could also been seen as the 
	// the current angle of the rover in the y-axis direction.
	float roverAngle = 0.0f;

	// the velocity array to use in intermediate calculations
	float roverVel[4] = { 0.0f, 0.0f, 0.0f, 0.0f };

	float minVel = 0.01f;
	float maxThetaVel = 0.6981f; // 2pi/9 rad/s == 40 degrees/s

	float rearCameraConfig[12] = { 0.0f };

	void setupMesh(MyMesh* mesh);
	void sendMatrices() const;

};

