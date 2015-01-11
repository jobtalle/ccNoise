#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>

#include <ccNoise/ccnoise.h>
#include <ccRandom/ccRandom.h>
#include <ccTrigonometry/ccTrigonometry.h>
#include <ccSort/ccSort.h>

#define _CCN_PERLIN_NORMALIZER 0.707107

typedef struct {
	int x, y;
} ccnPoint;

static float ccnInterpolate(float a, float b, float x, ccnInterpolationMethod interpolationMethod)
{
	switch(interpolationMethod) {
	case CCN_INTERP_LINEAR:
		return ccTriInterpolateLinear(a, b, x);
	case CCN_INTERP_QUADRATIC:
		return ccTriInterpolateQuadratic(a, b, x);
	case CCN_INTERP_QUADRATIC_INVERSE:
		return ccTriInterpolateQuadraticInverse(a, b, x);
	case CCN_INTERP_COSINE:
		return ccTriInterpolateCosine(a, b, x);
	case CCN_INTERP_PERLIN:
		return a + x*x*x*(x*(x * 6 - 15) + 10)*(b - a);
		break;
	default:
		return 0;
	}
}

static unsigned int ccnDistance(ccnPoint a, ccnPoint b, ccnDistanceMethod distanceMethod)
{
	switch(distanceMethod) {
	case CCN_DIST_MANHATTAN:
		return abs(a.x - b.x) + abs(a.y - b.y);
	case CCN_DIST_EUCLIDEAN:
		return (unsigned int)ccTriDistance(a.x, a.y, b.x, b.y);
	case CCN_DIST_CHEBYCHEV:
		return max(abs(a.x - b.x), abs(a.y - b.y));
	default:
		return 0;
	}
}

static int ccnWrapCoordinate(int coordinate, unsigned int period) {
	int positiveDeviation, negativeDeviation;

	if(period == CCN_INFINITE) return coordinate;

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

static void ccnStore(float *buffer, ccnStoreMethod method, float value)
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

int ccnGenerateWhiteNoise2D(
	ccnNoise *noise,
	ccnNoiseConfiguration *configuration)
{
	unsigned int size = noise->width * noise->height;
	unsigned int i;
	float multiplier = configuration->range.high - configuration->range.low;

	for(i = 0; i < size; i++) {
		int Y = i / noise->width;
		int X = i - Y * noise->width;
		ccnStore(noise->values + i, configuration->storeMethod, ccrGenerateFloatCoordinate(configuration->seed, X, Y) * multiplier + configuration->range.low);
	}

	return CCN_ERROR_NONE;
}

int ccnGenerateValueNoise2D(
	ccnNoise *noise,
	ccnNoiseConfiguration *configuration,
	unsigned int scale,
	ccnInterpolationMethod interpolationMethod)
{
	unsigned int size = noise->width * noise->height;
	unsigned int octaveWidth = noise->width / scale;
	unsigned int octaveHeight = noise->height / scale;
	unsigned int offsetHeight;
	unsigned int i, j, k;
	unsigned int yOffset = interpolationMethod == CCN_INTERP_CUBIC?1:0;

	float multiplier = configuration->range.high - configuration->range.low;
	float *xValues;

	ccnPoint offset;

	offsetHeight = octaveHeight + yOffset + (interpolationMethod == CCN_INTERP_CUBIC?2:1);
	xValues = malloc(noise->width * offsetHeight * sizeof(float));

	offset.x = configuration->tileConfiguration.xPeriod * octaveWidth;
	offset.y = configuration->tileConfiguration.yPeriod * octaveHeight;

	for(i = 0; i < offsetHeight; i++) {
		for(j = 0; j < noise->width; j++) {
			unsigned int octX = j / scale;

			float factor = (float)(j - octX * scale) / scale;

			if(interpolationMethod == CCN_INTERP_CUBIC) {
				float bufferedValues[4];

				if(factor == 0) {
					if(j == 0) {
						for(k = 0; k < 4; k++) {
							bufferedValues[k] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(octX - 1 + k + configuration->x * octaveWidth, offset.x), ccnWrapCoordinate(i + configuration->y * octaveHeight, offset.y));
						}
					}
					else {
						for(k = 0; k < 3; k++) {
							bufferedValues[k] = bufferedValues[k + 1];
						}

						bufferedValues[3] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(octX + 2 + configuration->x * octaveWidth, offset.x), ccnWrapCoordinate(i + configuration->y * octaveHeight, offset.y));
					}

					xValues[i * noise->width + j] = bufferedValues[1];
				}
				else {
					xValues[i * noise->width + j] = ccTriInterpolateCubic(bufferedValues[0], bufferedValues[1], bufferedValues[2], bufferedValues[3], factor);
				}
			}
			else {
				float bufferedValues[2];

				if(factor == 0) {
					if(j == 0) {
						bufferedValues[0] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(octX + configuration->x * octaveWidth, offset.x), ccnWrapCoordinate(i + configuration->y * octaveHeight, offset.y));
					}
					else {
						bufferedValues[0] = bufferedValues[1];
					}

					bufferedValues[1] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(octX + 1 + configuration->x * octaveWidth, offset.x), ccnWrapCoordinate(i + configuration->y * octaveHeight, offset.y));

					xValues[i * noise->width + j] = bufferedValues[0];
				}
				else {
					xValues[i * noise->width + j] = ccnInterpolate(bufferedValues[0], bufferedValues[1], factor, interpolationMethod);
				}
			}
		}
	}

	for(i = 0; i < size; i++) {
		unsigned int Y = i / noise->width;
		unsigned int octY = Y / scale;
		unsigned int index = (i - Y * noise->width) + (octY + yOffset) * noise->width;

		float factor = (float)(Y - octY * scale) / scale;

		if(factor == 0) {
			ccnStore(noise->values + i, configuration->storeMethod, xValues[index] * multiplier + configuration->range.low);
		}
		else {
			if(interpolationMethod == CCN_INTERP_CUBIC) {
				ccnStore(noise->values + i, configuration->storeMethod, ccTriInterpolateCubic(xValues[index - noise->width], xValues[index], xValues[index + noise->width], xValues[index + (noise->width << 1)], factor) * multiplier + configuration->range.low);
			}
			else {
				ccnStore(noise->values + i, configuration->storeMethod, ccnInterpolate(xValues[index], xValues[index + noise->width], factor, interpolationMethod) * multiplier + configuration->range.low);
			}
		}
	}
		
	free(xValues);

	return CCN_ERROR_NONE;
}

int ccnGenerateWorleyNoise2D(
	ccnNoise *noise,
	ccnNoiseConfiguration *configuration,
	unsigned int points,
	unsigned int n,
	int low, int high,
	ccnDistanceMethod distanceMethod,
	ccnInterpolationMethod interpolationMethod)
{
	unsigned int size = noise->width * noise->height;
	unsigned int i, j;
	unsigned int pointId = 0;
	unsigned int pointListSize = points * 9;
	ccnPoint offset;

	ccnPoint *pointList = malloc(pointListSize*sizeof(ccnPoint));
	int *pointsDistances = malloc(pointListSize*sizeof(unsigned int));

	unsigned int maxManhattanDistance = (unsigned int)(high * (2 / sqrt(2)));

	if(interpolationMethod == CCN_INTERP_CUBIC) return CCN_ERROR_INVALID_METHOD;

	if(configuration->tileConfiguration.tileMethod = CCN_TILE_NOT) configuration->tileConfiguration.xPeriod = configuration->tileConfiguration.yPeriod = CCN_INFINITE;

	for(offset.x = -1; offset.x <= 1; offset.x++) {
		for(offset.y = -1; offset.y <= 1; offset.y++) {
			ccRandomizer32 randomizer;

			ccrSeed32(&randomizer, ccrGenerateUintCoordinate(configuration->seed, ccnWrapCoordinate(configuration->x + offset.x, configuration->tileConfiguration.yPeriod), ccnWrapCoordinate(configuration->y + offset.y, configuration->tileConfiguration.yPeriod)));

			for(i = 0; i < points; i++) {
				pointList[pointId].x = (int)((ccrGenerateFloat32(&randomizer) + offset.x) * noise->width);
				pointList[pointId].y = (int)((ccrGenerateFloat32(&randomizer) + offset.y) * noise->height);
				pointId++;
			}
		}
	}

	for(i = 0; i < size; i++) {
		ccnPoint p;
		p.y = i / noise->width;
		p.x = i - p.y * noise->width;

		pointId = 0;
		for(j = 0; j < pointListSize; j++) {
			unsigned int manhattanDistance = ccnDistance(p, pointList[j], CCN_DIST_MANHATTAN);

			if(manhattanDistance < maxManhattanDistance) {
				if(distanceMethod == CCN_DIST_MANHATTAN) {
					pointsDistances[pointId] = manhattanDistance;
				}
				else {
					pointsDistances[pointId] = ccnDistance(p, pointList[j], distanceMethod);
				}

				pointId++;
			}
		}

		if(pointId > 1) ccsQuicksort(pointsDistances, 0, pointId);

		if(pointId <= n || pointsDistances[n] > high) {
			ccnStore(noise->values + i, configuration->storeMethod, configuration->range.high);
		}
		else if(pointsDistances[n] < low) {
			ccnStore(noise->values + i, configuration->storeMethod, configuration->range.low);
		}
		else {
			ccnStore(noise->values + i, configuration->storeMethod, ccnInterpolate(configuration->range.low, configuration->range.high, (float)(pointsDistances[n] - low) / (high - low), interpolationMethod));
		}
	}

	free(pointList);
	free(pointsDistances);

	return CCN_ERROR_NONE;
}

int ccnGeneratePerlinNoise2D(
	ccnNoise *noise,
	ccnNoiseConfiguration *configuration,
	unsigned int scale,
	ccnInterpolationMethod interpolationMethod)
{
	unsigned int xSteps = (noise->width / scale) + 1;
	unsigned int ySteps = noise->height / scale;
	unsigned int totalSteps = xSteps * (ySteps + 1);
	unsigned int size = noise->width * noise->height;
	unsigned int i;

	float multiplier = configuration->range.high - configuration->range.low;
	float *vectors = malloc(sizeof(float)* (totalSteps << 1));
	
	ccnPoint offset = (ccnPoint){ configuration->x * (xSteps - 1), configuration->y * ySteps };

	if(configuration->tileConfiguration.tileMethod = CCN_TILE_NOT) {
		configuration->tileConfiguration.xPeriod = configuration->tileConfiguration.yPeriod = CCN_INFINITE;
	}
	else{
		configuration->tileConfiguration.xPeriod *= xSteps - 1;
		configuration->tileConfiguration.yPeriod *= ySteps;
	}

	for(i = 0; i < totalSteps; i++) {
		int Y = i / xSteps;
		float radians = (float)(ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(i - Y * xSteps + offset.x, configuration->tileConfiguration.xPeriod), ccnWrapCoordinate(Y + offset.y, configuration->tileConfiguration.yPeriod)) * CC_TRI_PI_DOUBLE);

		vectors[i << 1] = (float)cos(radians);
		vectors[(i << 1) + 1] = (float)sin(radians);
	}

	for(i = 0; i < size; i++) {
		int Y = i / noise->width;
		int X = i - Y * noise->width;
		int xStep = X / scale;
		int yStep = Y / scale;

		unsigned int indexTop = (xStep + yStep * xSteps) << 1;
		unsigned int indexBottom = indexTop + (xSteps << 1);

		float factorX = (float)(X - xStep * scale) / scale;

		float vecX = (float)(X - xStep * scale) / scale;
		float vecY = (float)(Y - yStep * scale) / scale;

		ccnStore(noise->values + i, configuration->storeMethod,
			(float)((ccnInterpolate(
			ccnInterpolate(
			vectors[indexTop] * vecX + vectors[indexTop + 1] * vecY,
			vectors[indexTop + 2] * (vecX - 1.0f) + vectors[indexTop + 3] * vecY,
			factorX, interpolationMethod),
			ccnInterpolate(
			vectors[indexBottom] * vecX + vectors[indexBottom + 1] * (vecY - 1.0f),
			vectors[indexBottom + 2] * (vecX - 1.0f) + vectors[indexBottom + 3] * (vecY - 1.0f),
			factorX, interpolationMethod),
			(float)(Y - yStep * scale) / scale, interpolationMethod) + _CCN_PERLIN_NORMALIZER) * _CCN_PERLIN_NORMALIZER * multiplier) + configuration->range.low);
	}

	free(vectors);

	return CCN_ERROR_NONE;
}