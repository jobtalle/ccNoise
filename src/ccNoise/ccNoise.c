#include <stdlib.h>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>
#include <ccTrigonometry/ccTrigonometry.h>
#include <ccSort/ccSort.h>

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
	default:
		return 0;
	}
}

static unsigned int ccnDistance(ccPoint a, ccPoint b, ccnDistanceMethod distanceMethod)
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

int ccnGenerateWorleyNoise(
	float **buffer,
	unsigned int seed,
	ccnTileConfiguration *tileConfig,
	int x, int y,
	unsigned int width, unsigned int height,
	ccnStoreMethod storeMethod,
	ccnRange range,

	unsigned int points,
	unsigned int n,
	int low, int high,
	ccnDistanceMethod distanceMethod,
	ccnInterpolationMethod interpolationMethod)
{
	unsigned int size = width * height;
	unsigned int i, j;
	unsigned int pointId = 0;
	unsigned int pointListSize = points * 9;
	ccPoint offset;

	ccPoint *pointList = malloc(pointListSize*sizeof(ccPoint));
	int *pointsDistances = malloc(pointListSize*sizeof(unsigned int));

	unsigned int maxManhattanDistance = (unsigned int)(high * (2 / sqrt(2)));

	if(interpolationMethod == CCN_INTERP_CUBIC) return CCN_ERROR_INVALID_METHOD;

	if(tileConfig->tileMethod = CCN_TILE_NOT) tileConfig->xPeriod = tileConfig->yPeriod = CCN_INFINITE;

	for(offset.x = -1; offset.x <= 1; offset.x++) {
		for(offset.y = -1; offset.y <= 1; offset.y++) {
			ccRandomizer32 randomizer;

			ccrSeed32(&randomizer, ccrGenerateUintCoordinate(seed, ccnWrapCoordinate(x + offset.x, tileConfig->xPeriod), ccnWrapCoordinate(y + offset.y, tileConfig->yPeriod)));

			for(i = 0; i < points; i++) {
				pointList[pointId].x = (int)((ccrGenerateFloat32(&randomizer) + offset.x) * width);
				pointList[pointId].y = (int)((ccrGenerateFloat32(&randomizer) + offset.y) * height);
				pointId++;
			}
		}
	}

	for(i = 0; i < size; i++) {
		ccPoint p;
		p.y = i / width;
		p.x = i - p.y * width;

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
			ccnStore(*buffer + i, storeMethod, range.high);
		}
		else if(pointsDistances[n] < low) {
			ccnStore(*buffer + i, storeMethod, range.low);
		}
		else {
			ccnStore(*buffer + i, storeMethod, ccnInterpolate(range.low, range.high, (float)(pointsDistances[n] - low) / (high - low), interpolationMethod));
		}
	}

	free(pointList);
	free(pointsDistances);

	return CCN_ERROR_NONE;
}

int ccnGenerateWhiteNoise(
	float **buffer,
	unsigned int seed,
	unsigned int width, unsigned int height,
	ccnStoreMethod storeMethod,
	ccnRange range)
{
	unsigned int size = width * height;
	unsigned int i;
	float multiplier = range.high - range.low;
	
	ccRandomizer32 randomizer;

	ccrSeed32(&randomizer, seed);

	for(i = 0; i < size; i++) {
		ccnStore(*buffer + i, storeMethod, ccrGenerateFloat32(&randomizer) * multiplier + range.low);
	}

	return CCN_ERROR_NONE;
}

int ccnGenerateValueNoise(
	float **buffer,
	unsigned int seed,
	ccnTileConfiguration *tileConfig,
	int x, int y,
	unsigned int width, unsigned int height,
	ccnStoreMethod storeMethod,
	ccnRange range,

	unsigned int scale,
	ccnInterpolationMethod interpolationMethod)
{
	unsigned int size = width * height;
	unsigned int octaveWidth = width / scale;
	unsigned int octaveHeight = height / scale;
	unsigned int offsetHeight;
	unsigned int i, j, k;
	unsigned int yOffset = interpolationMethod == CCN_INTERP_CUBIC?1:0;

	float multiplier = range.high - range.low;
	float *xValues;

	ccPoint offset;

	offsetHeight = octaveHeight + yOffset + (interpolationMethod == CCN_INTERP_CUBIC?2:1);
	xValues = malloc(width * offsetHeight * sizeof(float));

	offset.x = tileConfig->xPeriod * octaveWidth;
	offset.y = tileConfig->yPeriod * octaveHeight;

	for(i = 0; i < offsetHeight; i++) {
		for(j = 0; j < width; j++) {
			unsigned int octX = j / scale;

			float factor = (float)(j - octX * scale) / scale;

			if(interpolationMethod == CCN_INTERP_CUBIC) {
				float bufferedValues[4];

				if(factor == 0) {
					if(j == 0) {
						for(k = 0; k < 4; k++) {
							bufferedValues[k] = ccrGenerateFloatCoordinate(seed, ccnWrapCoordinate(octX - 1 + k + x * octaveWidth, offset.x), ccnWrapCoordinate(i + y * octaveHeight, offset.y));
						}
					}
					else {
						for(k = 0; k < 3; k++) {
							bufferedValues[k] = bufferedValues[k + 1];
						}

						bufferedValues[3] = ccrGenerateFloatCoordinate(seed, ccnWrapCoordinate(octX + 2 + x * octaveWidth, offset.x), ccnWrapCoordinate(i + y * octaveHeight, offset.y));
					}

					xValues[i * width + j] = bufferedValues[1];
				}
				else {
					xValues[i * width + j] = ccTriInterpolateCubic(bufferedValues[0], bufferedValues[1], bufferedValues[2], bufferedValues[3], factor);
				}
			}
			else {
				float bufferedValues[2];

				if(factor == 0) {
					if(j == 0) {
						bufferedValues[0] = ccrGenerateFloatCoordinate(seed, ccnWrapCoordinate(octX + x * octaveWidth, offset.x), ccnWrapCoordinate(i + y * octaveHeight, offset.y));
					}
					else {
						bufferedValues[0] = bufferedValues[1];
					}

					bufferedValues[1] = ccrGenerateFloatCoordinate(seed, ccnWrapCoordinate(octX + 1 + x * octaveWidth, offset.x), ccnWrapCoordinate(i + y * octaveHeight, offset.y));

					xValues[i * width + j] = bufferedValues[0];
				}
				else {
					xValues[i * width + j] = ccnInterpolate(bufferedValues[0], bufferedValues[1], factor, interpolationMethod);
				}
			}
		}
	}

	for(i = 0; i < size; i++) {
		unsigned int Y = i / width;
		unsigned int octY = Y / scale;
		unsigned int index = (i - Y * width) + (octY + yOffset) * width;

		float factor = (float)(Y - octY * scale) / scale;

		if(factor == 0) {
			ccnStore(*buffer + i, storeMethod, xValues[index] * multiplier + range.low);
		}
		else {
			if(interpolationMethod == CCN_INTERP_CUBIC) {
				ccnStore(*buffer + i, storeMethod, ccTriInterpolateCubic(xValues[index - width], xValues[index], xValues[index + width], xValues[index + (width << 1)], factor) * multiplier + range.low);
			}
			else {
				ccnStore(*buffer + i, storeMethod, ccnInterpolate(xValues[index], xValues[index + width], factor, interpolationMethod) * multiplier + range.low);
			}
		}
	}
		
	free(xValues);

	return CCN_ERROR_NONE;
}