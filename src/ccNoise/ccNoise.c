#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>

#include <ccNoise/ccnoise.h>
#include <ccRandom/ccRandom.h>
#include <ccTrigonometry/ccTrigonometry.h>
#include <ccSort/ccSort.h>

#define _CCN_PERLIN_NORMALIZER         0.707107
#define _CCN_MANHATTAN_DISTANCE_FACTOR 1.414214

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

static int ccnFloorMod(int x, int y) {
	return x >= 0?x % y:y + ((x + 1) % y) - 1;
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
	unsigned int xSteps = (unsigned int)ceil((float)noise->width / scale);
	unsigned int ySteps = (unsigned int)ceil((float)noise->width / scale);
	unsigned int offsetHeight;
	unsigned int xPeriod;
	unsigned int yPeriod;
	unsigned int xOffset = 0;
	unsigned int yOffset = 0;
	unsigned int interpolationYoffset = interpolationMethod == CCN_INTERP_CUBIC?1:0;
	unsigned int i, j, k;

	float multiplier = configuration->range.high - configuration->range.low;
	float *xValues;

	ccnPoint offset = (ccnPoint){ configuration->x * xSteps, configuration->y * ySteps };

#ifdef _DEBUG
	if(scale & (scale - 1)) return CCN_ERROR_NO_POWER_OF_2;
#endif

	if(noise->width < scale) {
		offset.x = (int)floor(offset.x * ((float)noise->width / scale));
		xOffset = ccnFloorMod(configuration->x, scale / noise->width) * noise->width;
	}

	if(noise->height < scale) {
		offset.y = (int)floor(offset.y * ((float)noise->height / scale));
		yOffset = ccnFloorMod(configuration->y, scale / noise->height) * noise->height;
	}

	if(configuration->tileConfiguration.tileMethod = CCN_TILE_NOT) {
		xPeriod = yPeriod = CCN_INFINITE;
	}
	else{
		xPeriod = (unsigned int)(configuration->tileConfiguration.xPeriod * ((float)noise->width / scale));
		yPeriod = (unsigned int)(configuration->tileConfiguration.yPeriod * ((float)noise->height / scale));
	}

	offsetHeight = ySteps + (interpolationMethod == CCN_INTERP_CUBIC?3:1);
	xValues = malloc(noise->width * offsetHeight * sizeof(float));

	for(i = 0; i < offsetHeight; i++) {
		for(j = 0; j < noise->width; j++) {
			unsigned int octX = j / scale;

			float factor = (float)(j - octX * scale) / scale;
			float interpFactor = factor + (float)xOffset / scale;

			if(interpolationMethod == CCN_INTERP_CUBIC) {
				float bufferedValues[4];

				if(factor == 0) {
					if(j == 0) {
						for(k = 0; k < 4; k++) {
							bufferedValues[k] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(octX - 1 + k + offset.x, xPeriod), ccnWrapCoordinate(i + offset.y, yPeriod));
						}
					}
					else {
						for(k = 0; k < 3; k++) {
							bufferedValues[k] = bufferedValues[k + 1];
						}

						bufferedValues[3] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(octX + 2 + offset.x, xPeriod), ccnWrapCoordinate(i + offset.y, yPeriod));
					}
				}

				xValues[i * noise->width + j] = ccTriInterpolateCubic(bufferedValues[0], bufferedValues[1], bufferedValues[2], bufferedValues[3], interpFactor);
			}
			else {
				float bufferedValues[2];

				if(factor == 0) {
					if(j == 0) {
						bufferedValues[0] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(octX + offset.x, xPeriod), ccnWrapCoordinate(i + offset.y, yPeriod));
					}
					else {
						bufferedValues[0] = bufferedValues[1];
					}

					bufferedValues[1] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(octX + 1 + offset.x, xPeriod), ccnWrapCoordinate(i + offset.y, yPeriod));
				}
				
				xValues[i * noise->width + j] = ccnInterpolate(bufferedValues[0], bufferedValues[1], interpFactor, interpolationMethod);
			}
		}
	}

	for(i = 0; i < size; i++) {
		unsigned int Y = i / noise->width;
		unsigned int octY = Y / scale;
		unsigned int index = (i - Y * noise->width) + (octY + interpolationYoffset) * noise->width;

		float factor = (float)(Y + yOffset - octY * scale) / scale;

		if(interpolationMethod == CCN_INTERP_CUBIC) {
			if(factor == 0) {
				ccnStore(noise->values + i, configuration->storeMethod, (xValues[index] * .5f + 0.25f) * multiplier + configuration->range.low);
			}
			else {
				ccnStore(noise->values + i, configuration->storeMethod, (ccTriInterpolateCubic(xValues[index - noise->width], xValues[index], xValues[index + noise->width], xValues[index + (noise->width << 1)], factor) * .5f + 0.25f) * multiplier + configuration->range.low);
			}
		}
		else {
			if(factor == 0) {
				ccnStore(noise->values + i, configuration->storeMethod, xValues[index] * multiplier + configuration->range.low);
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
	unsigned int xPeriod = configuration->tileConfiguration.xPeriod;
	unsigned int yPeriod = configuration->tileConfiguration.yPeriod;
	ccnPoint offset;

	ccnPoint *pointList = malloc(pointListSize*sizeof(ccnPoint));
	int *pointsDistances = malloc(pointListSize*sizeof(unsigned int));

	unsigned int maxManhattanDistance = (unsigned int)(high * _CCN_MANHATTAN_DISTANCE_FACTOR);

#ifdef _DEBUG
	if(interpolationMethod == CCN_INTERP_CUBIC) return CCN_ERROR_INVALID_METHOD;
#endif

	if(configuration->tileConfiguration.tileMethod = CCN_TILE_NOT) xPeriod = yPeriod = CCN_INFINITE;

	for(offset.x = -1; offset.x <= 1; offset.x++) {
		for(offset.y = -1; offset.y <= 1; offset.y++) {
			ccRandomizer32 randomizer;

			ccrSeed32(&randomizer, ccrGenerateUintCoordinate(configuration->seed, ccnWrapCoordinate(configuration->x + offset.x, xPeriod), ccnWrapCoordinate(configuration->y + offset.y, yPeriod)));

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
	unsigned int xSteps = (unsigned int)ceil((float)noise->width / scale);
	unsigned int ySteps = (unsigned int)ceil((float)noise->height / scale);
	unsigned int totalSteps = (xSteps + 1) * (ySteps + 1);
	unsigned int size = noise->width * noise->height;
	unsigned int xPeriod;
	unsigned int yPeriod;
	unsigned int xOffset = 0;
	unsigned int yOffset = 0;
	unsigned int i;

	float multiplier = configuration->range.high - configuration->range.low;
	float *vectors;

	ccnPoint offset = (ccnPoint){ configuration->x * xSteps, configuration->y * ySteps };

#ifdef _DEBUG
	if(interpolationMethod == CCN_INTERP_CUBIC) return CCN_ERROR_INVALID_METHOD;
	if(scale & (scale - 1)) return CCN_ERROR_NO_POWER_OF_2;
#endif

	if(noise->width < scale) {
		offset.x = (int)floor(offset.x * ((float)noise->width / scale));
		xOffset = ccnFloorMod(configuration->x, scale / noise->width) * noise->width;
	}

	if(noise->height < scale) {
		offset.y = (int)floor(offset.y * ((float)noise->height / scale));
		yOffset = ccnFloorMod(configuration->y, scale / noise->height) * noise->height;
	}

	vectors = malloc(sizeof(float)* (totalSteps << 1));

	if(configuration->tileConfiguration.tileMethod = CCN_TILE_NOT) {
		xPeriod = yPeriod = CCN_INFINITE;
	}
	else{
		xPeriod = (unsigned int)(configuration->tileConfiguration.xPeriod * ((float)noise->width / scale));
		yPeriod = (unsigned int)(configuration->tileConfiguration.yPeriod * ((float)noise->height / scale));
	}

	for(i = 0; i < totalSteps; i++) {
		int Y = i / (xSteps + 1);
		int X = i - Y * (xSteps + 1);
		float radians = (float)(ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(X + offset.x, xPeriod), ccnWrapCoordinate(Y + offset.y, yPeriod)) * CC_TRI_PI_DOUBLE);

		vectors[i << 1] = (float)cos(radians);
		vectors[(i << 1) + 1] = (float)sin(radians);
	}

	for(i = 0; i < size; i++) {
		int Y = i / noise->width;
		int X = i - Y * noise->width;
		int xStep = X / scale;
		int yStep = Y / scale;

		X += xOffset;
		Y += yOffset;

		unsigned int indexTop = (xStep + yStep * (xSteps + 1)) << 1;
		unsigned int indexBottom = indexTop + ((xSteps + 1) << 1);

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