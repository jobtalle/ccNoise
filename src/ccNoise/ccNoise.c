#include <math.h>

#include <ccNoise/ccNoise.h>
#include <ccTrigonometry/ccTrigonometry.h>

#include "ccNoiseInternal.h"

#ifndef max
#define max(x, y) (x)>(y)?(x):(y)
#endif

#ifdef __GNUC__
#define max(a,b) \
	({ __typeof__ (a) _a = (a); \
	__typeof__ (b) _b = (b); \
	_a > _b ? _a : _b; })
#endif

float ccnInterpolate(float a, float b, float x, ccnInterpolationMethod interpolationMethod)
{
	switch(interpolationMethod) {
	case CCN_INTERP_LINEAR:
		return a * (1 - x) + b*x;
	case CCN_INTERP_COSINE:
	{
		float factor = (1.0f - cosf(x * (float)CC_TRI_PI)) * .5f;
		return a * (1.0f - factor) + b * factor;
	}
	case CCN_INTERP_PERLIN:
		return a + x*x*x*(x*(x * 6 - 15) + 10)*(b - a);
		break;
	default:
		return 0;
	}
}

float ccnInterpolateCubic(float a, float b, float c, float d, float x)
{
	float p = (d - c) - (a - b);
	return ccTriCubed(x) * p + ccTriSquared(x) * ((a - b) - p) + x * (c - a) + b;
}

unsigned int ccnDistance(ccnPoint a, ccnPoint b, ccnDistanceMethod distanceMethod)
{
	switch(distanceMethod) {
	case CCN_DIST_MANHATTAN:
		return abs(a.x - b.x) + abs(a.y - b.y);
	case CCN_DIST_EUCLIDEAN:
		return (unsigned int)ccTriDistance((float)a.x, (float)a.y, (float)b.x, (float)b.y);
	case CCN_DIST_CHEBYCHEV:
		return max(abs(a.x - b.x), abs(a.y - b.y));
	default:
		return 0;
	}
}

int ccnWrapCoordinate(int coordinate, unsigned int period) {
	int positiveDeviation, negativeDeviation;

	if(period == CCN_INFINITE || period == 0) return coordinate;

	positiveDeviation = period >> 1;
	negativeDeviation = -(int)(period - positiveDeviation);

	if(coordinate > positiveDeviation) {
		coordinate -= ((coordinate - negativeDeviation - 1) / period) * period;
	}
	else if(coordinate <= negativeDeviation) {
		coordinate += ((-coordinate + positiveDeviation) / period) * period;
	}

	return coordinate;
}

void ccnStore(float *buffer, ccnStoreMethod method, float value)
{
	switch(method) {
	case CCN_STORE_SET:
		*buffer = value;
		break;
	case CCN_STORE_ADD:
		*buffer += value;
		break;
	case CCN_STORE_SUBTRACT:
		*buffer -= value;
		break;
	case CCN_STORE_MULTIPLY:
		*buffer *= value;
		break;
	case CCN_STORE_DIVIDE:
		*buffer /= value;
		break;
	}
}

int ccnFloorMod(int x, int y) {
	return x >= 0?x % y:y + ((x + 1) % y) - 1;
}