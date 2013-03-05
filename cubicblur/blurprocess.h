
#ifndef BLUR_PROCESS_H
#define BLUR_PROCESS_H

typedef struct cubepixel {
	float r;
	float g;
	float b;

	float x;
	float y;
	float z;
} CubePixel;

typedef struct pixel {
	float r;
	float g;
	float b;
} Pixel;

float** blurFaces( float ** faces, int inputSize, int outputSize, float power, float curve );

float** computeBlur( CubePixel ** faces, int outputSize, float power, float curve );

CubePixel** createCpMap( float ** faces, int inputSize );

#endif BLUR_PROCESS_H