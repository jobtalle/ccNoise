#pragma once

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>

typedef struct {
	int x, y;
} ccnPoint;

float ccnInterpolate(float a, float b, float x, ccnInterpolationMethod interpolationMethod);
float ccnInterpolateCubic(float a, float b, float c, float d, float x);
unsigned int ccnDistance(ccnPoint a, ccnPoint b, ccnDistanceMethod distanceMethod);
int ccnWrapCoordinate(int coordinate, unsigned int period);
void ccnStore(float *buffer, ccnStoreMethod method, float value);
int ccnFloorMod(int x, int y);