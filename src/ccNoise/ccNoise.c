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

int ccnGenerateWorleyNoise(
	float **buffer,
	unsigned int seed,
	ccnTileConfiguration *tileConfig,
	int x, int y,
	unsigned int width, unsigned int height,
	unsigned int points,
	unsigned int n,
	int low, int high,
	float lowValue, float highValue,
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
		
		if(pointId <= n) {
			(*buffer)[i] = highValue;
		}
		else if(pointsDistances[n] > high) {
			(*buffer)[i] = highValue;
		}
		else if(pointsDistances[n] < low) {
			(*buffer)[i] = lowValue;
		}
		else {
			(*buffer)[i] = ccnInterpolate(lowValue, highValue, (float)(pointsDistances[n] - low) / (high - low), interpolationMethod);
		}
	}

	free(pointList);
	free(pointsDistances);

	return CCN_ERROR_NONE;
}

int ccnGenerateWhiteNoise(
	float **buffer,
	unsigned int seed,
	unsigned int width, unsigned int height)
{
	unsigned int size = width * height;
	unsigned int i;
	
	ccRandomizer32 randomizer;

	ccrSeed32(&randomizer, seed);

	for(i = 0; i < size; i++) {
		(*buffer)[i] = ccrGenerateFloat32(&randomizer);
	}

	return CCN_ERROR_NONE;
}

int ccnGenerateValueNoise(
	float **buffer,
	unsigned int seed,
	ccnTileConfiguration *tileConfig,
	int x, int y,
	unsigned int width, unsigned int height,
	unsigned int octaves,
	unsigned int maxOctave,
	ccnInterpolationMethod interpolationMethod)
{
	unsigned int size = width * height;
	unsigned int octaveSize = maxOctave;
	unsigned int i, j, k, l;
	
	float influence = 0.5f;

	ccPoint negativeOffset;
	ccPoint positiveOffset;

	// TODO: inline the following
	if(interpolationMethod == CCN_INTERP_CUBIC) {
		negativeOffset.x = negativeOffset.y = 1;
		positiveOffset.x = positiveOffset.y = 2;
	}
	else {
		negativeOffset.x = negativeOffset.y = 0;
		positiveOffset.x = positiveOffset.y = 1;
	}

	octaveSize = maxOctave;

	for(i = 0; i < size; i++) {
		(*buffer)[i] = 0;
	}

	for(i = 0; i < octaves; i++) {
		unsigned int octaveWidth = width / octaveSize;
		unsigned int octaveHeight = height / octaveSize;
		unsigned int offsetHeight = octaveHeight + negativeOffset.y + positiveOffset.y;
		ccPoint offset;

		float *xValues = malloc(width * offsetHeight * sizeof(float));

		offset.x = tileConfig->xPeriod * octaveWidth;
		offset.y = tileConfig->yPeriod * octaveHeight;

		for(j = 0; j < offsetHeight; j++) {
			for(k = 0; k < width; k++) {
				unsigned int octX = k / octaveSize;

				float factor = (float)(k - octX * octaveSize) / octaveSize;

				if(interpolationMethod == CCN_INTERP_CUBIC) {
					float bufferedValues[4];

					if(factor == 0) {
						if(k == 0) {
							for(l = 0; l < 4; l++) {
								bufferedValues[l] = ccrGenerateFloatCoordinate(seed, ccnWrapCoordinate(octX - 1 + l + x * octaveWidth, offset.x), ccnWrapCoordinate(j + y * octaveHeight, offset.y));
							}
						}
						else {
							for(l = 0; l < 3; l++) {
								bufferedValues[l] = bufferedValues[l + 1];
							}

							bufferedValues[3] = ccrGenerateFloatCoordinate(seed, ccnWrapCoordinate(octX + 2 + x * octaveWidth, offset.x), ccnWrapCoordinate(j + y * octaveHeight, offset.y));
						}

						xValues[j * width + k] = bufferedValues[1];
					}
					else {
						xValues[j * width + k] = ccTriInterpolateCubic(bufferedValues[0], bufferedValues[1], bufferedValues[2], bufferedValues[3], factor);
					}
				}
				else {
					float bufferedValues[2];

					if(factor == 0) {
						if(k == 0) {
							bufferedValues[0] = ccrGenerateFloatCoordinate(seed, ccnWrapCoordinate(octX + x * octaveWidth, offset.x), ccnWrapCoordinate(j + y * octaveHeight, offset.y));
						}
						else {
							bufferedValues[0] = bufferedValues[1];
						}

						bufferedValues[1] = ccrGenerateFloatCoordinate(seed, ccnWrapCoordinate(octX + 1 + x * octaveWidth, offset.x), ccnWrapCoordinate(j + y * octaveHeight, offset.y));

						xValues[j * width + k] = bufferedValues[0];
					}
					else {
						xValues[j * width + k] = ccnInterpolate(bufferedValues[0], bufferedValues[1], factor, interpolationMethod);
					}
				}
			}
		}

		for(j = 0; j < size; j++) {
			unsigned int Y = j / width;
			unsigned int octY = Y / octaveSize;
			unsigned int index = (j - Y * width) + (octY + negativeOffset.y) * width;

			float factor = (float)(Y - octY * octaveSize) / octaveSize;

			if(factor == 0) {
				(*buffer)[j] += xValues[index] * influence;
			}
			else {
				if(interpolationMethod == CCN_INTERP_CUBIC) {
					(*buffer)[j] += ccTriInterpolateCubic(xValues[index - width], xValues[index], xValues[index + width], xValues[index + (width << 1)], factor) * influence;
				}
				else {
					(*buffer)[j] += ccnInterpolate(xValues[index], xValues[index + width], factor, interpolationMethod) * influence;
				}
			}
		}
		
		free(xValues);

		influence *= 0.5f;

		octaveSize >>= 1;
		seed++;
		if(octaveSize == 0) {
			if(octaves == CCN_INFINITE) {
				break;
			}
			else {
				return CCN_ERROR_INVALID_ARGUMENT_RANGE;
			}
		}
	}

	return CCN_ERROR_NONE;
}