#pragma once

#define CAMERA_FIXED_ORTHOGNAL 0
#define CAMERA_FIXED_PERSPECTIVE 1
#define CAMERA_FOLLOWING 2

class Camera {
public:
	int mode = 0;
	float camPos[3] = { 0.0f, 0.0f, 0.0f };
	float camTarget[3] = { 0.0f, 0.0f, 0.0f };
};