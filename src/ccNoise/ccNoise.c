#include <stdlib.h>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>
#include <ccTrigonometry/ccTrigonometry.h>
#include <ccSort/ccSort.h>

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

static void ccnGenerateOffsetNoise(
	float **buffer,
	unsigned int seed,
	int x, int y,
	unsigned int width, unsigned int height,
	ccnPoint negativeOffset, ccnPoint positiveOffset)
{
	unsigned int X, Y, offset;

	unsigned int totalWidth = width + negativeOffset.x + positiveOffset.x;
	unsigned int totalHeight = height + negativeOffset.y + positiveOffset.y;
	unsigned int totalSize = totalWidth * totalHeight;

	float *whiteNoiseBuffer;

	*buffer = malloc(sizeof(float)*totalSize);

	// Generate central noise
	ccnGenerateWhiteNoise(&whiteNoiseBuffer, seed, width, height);

	for(Y = 0; Y < height; Y++) {
		for(X = 0; X < width; X++) {
			(*buffer)[(Y + negativeOffset.y) * totalWidth + X + negativeOffset.x] = whiteNoiseBuffer[(Y * width) + X];
		}
	}

	free(whiteNoiseBuffer);
}

unsigned int ccnCoordinateUid(int x, int y)
{
	unsigned int shell = max(abs(x), abs(y));
	unsigned int uid = (x + y + (shell << 1)) << 1;
	
	if(shell != 0) {
		uid += ccTriSquared((shell << 1) - 1);
		if(x > y || (x == shell && y == shell)) {
			uid--;
		}
	}

	return uid;
}

int ccnGenerateWorleyNoise(
	float **buffer,
	unsigned int seed,
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
	ccnPoint offset;

	ccnPoint *pointList = malloc(pointListSize*sizeof(ccnPoint));
	int *pointsDistances = malloc(pointListSize*sizeof(unsigned int));

	unsigned int maxManhattanDistance = (unsigned int)(high * (2 / sqrt(2)));

	if(interpolationMethod == CCN_INTERP_CUBIC) return CCN_ERROR_INVALID_METHOD;

	*buffer = malloc(sizeof(float)*size);

	for(offset.x = -1; offset.x <= 1; offset.x++) {
		for(offset.y = -1; offset.y <= 1; offset.y++) {
			ccRandomizer32 randomizer;

			ccrSeed32(&randomizer, seed + ccnCoordinateUid(x + offset.x, y + offset.y));

			for(i = 0; i < points; i++) {
				pointList[pointId].x = (int)((ccrGenerateFloat32(&randomizer) + offset.x) * width);
				pointList[pointId].y = (int)((ccrGenerateFloat32(&randomizer) + offset.y) * height);
				pointId++;
			}
		}
	}

	for(i = 0; i < size; i++) {
		ccnPoint p;
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

	*buffer = malloc(size * sizeof(float));

	for(i = 0; i < size; i++) {
		(*buffer)[i] = ccrGenerateFloat32(&randomizer);
	}

	return CCN_ERROR_NONE;
}

int ccnGenerateValueNoise(
	float **buffer,
	unsigned int seed,
	ccnFlags tileFlags,
	int hPeriod, int vPeriod,
	int x, int y,
	unsigned int width, unsigned int height,
	unsigned int octaves,
	unsigned int maxOctave,
	ccnInterpolationMethod interpolationMethod)
{
	unsigned int size = width * height;
	unsigned int octaveSize = maxOctave;
	unsigned int i, j, k;
	
	float influence = 0.5f;

	ccnPoint negativeOffset;
	ccnPoint positiveOffset;

	if(interpolationMethod == CCN_INTERP_CUBIC) {
		negativeOffset.x = 1;
		negativeOffset.y = 1;
		positiveOffset.x = 2;
		positiveOffset.y = 2;
	}
	else {
		negativeOffset.x = 0;
		negativeOffset.y = 0;
		positiveOffset.x = 1;
		positiveOffset.y = 1;
	}

	octaveSize = maxOctave;

	*buffer = calloc(size, sizeof(float));

	for(i = 0; i < octaves; i++) {
		unsigned int xSteps = width / octaveSize;
		unsigned int ySteps = height / octaveSize;
		unsigned int offsetWidth = xSteps + negativeOffset.x + positiveOffset.x;
		unsigned int offsetHeight = ySteps + negativeOffset.y + positiveOffset.y;

		float *offsetNoise;
		float *xValues = malloc(width * offsetHeight * sizeof(float));

		ccnGenerateOffsetNoise(&offsetNoise, seed, x, y, xSteps, ySteps, negativeOffset, positiveOffset);

		// Find x values

		for(j = 0; j < offsetHeight; j++) {
			for(k = 0; k < width; k++) {
				unsigned int octX = k / octaveSize;
				unsigned int offsetIndex = j * offsetWidth + negativeOffset.x + octX;

				float factor = (float)(k - octX * octaveSize) / octaveSize;

				if(factor == 0) {
					xValues[j * width + k] = offsetNoise[offsetIndex];
				}
				else {
					if(interpolationMethod == CCN_INTERP_CUBIC) {
						xValues[j * width + k] = ccTriInterpolateCubic(offsetNoise[offsetIndex - 1], offsetNoise[offsetIndex], offsetNoise[offsetIndex + 1], offsetNoise[offsetIndex + 2], factor);
					}
					else {
						xValues[j * width + k] = ccnInterpolate(offsetNoise[offsetIndex], offsetNoise[offsetIndex + 1], factor, interpolationMethod);
					}
				}
			}
		}

		// Interpolate x values

		for(j = 0; j < size; j++) {
			unsigned int Y = j / width;
			unsigned int X = j - Y * width;

			unsigned octY = Y / octaveSize;

			float factor = (float)(Y - octY * octaveSize) / octaveSize;

			//TODO: handle factor==0 case

			unsigned int index = X + (octY + negativeOffset.y) * width;

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