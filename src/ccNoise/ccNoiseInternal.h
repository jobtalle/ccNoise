#pragma once

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>

typedef struct {
	int32_t x, y;
} ccnPoint;

float ccnInterpolate(float a, float b, float x, ccnInterpolationMethod interpolationMethod);
float ccnInterpolateCubic(float a, float b, float c, float d, float x);
uint32_t ccnDistance(ccnPoint a, ccnPoint b, ccnDistanceMethod distanceMethod);
int32_t ccnWrapCoordinate(int32_t coordinate, uint32_t period);
void ccnStore(float *buffer, ccnStoreMethod method, float value);
int32_t ccnFloorMod(int32_tx, int32_ty);
