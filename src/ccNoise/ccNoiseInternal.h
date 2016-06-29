#pragma once

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>

typedef struct {
	int32_t x, y;
} ccnPoint;

float ccnInterpolate(const float a, const float b, const float x, const ccnInterpolationMethod interpolationMethod);
float ccnInterpolateCubic(const float a, const float b, const float c, const float d, const float x);
uint32_t ccnDistance(const ccnPoint a, const ccnPoint b, const ccnDistanceMethod distanceMethod);
int32_t ccnWrapCoordinate(int32_t coordinate, const uint32_t period);
void ccnStore(float *buffer, const ccnStoreMethod method, const float value);
int32_t ccnFloorMod(const int32_t x, const int32_t y);
