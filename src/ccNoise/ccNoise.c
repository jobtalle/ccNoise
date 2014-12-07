#include <stdlib.h>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>
#include <ccTrigonometry/ccTrigonometry.h>
#include <ccSort/ccSort.h>

#define _CCN_CEIL_DIV_INT(n, div) ((n + div - 1) / div)

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

int ccnGenerateFractalNoise(
	float **buffer,
	unsigned int seed,
	bool makeTileable,
	int x, int y,
	unsigned int width, unsigned int height,
	unsigned int octaves,
	unsigned int maxOctave,
	ccnInterpolationMethod interpolationMethod)
{
	unsigned int size = width * height;
	unsigned int octaveSize = maxOctave;
	unsigned int i, j, k;
	
	float influence = 0.5f; // TODO: pick number to make range [0, 1]

	ccRandomizer32 randomizer;

	octaveSize = maxOctave;

	*buffer = calloc(size, sizeof(float));

	for(i = 0; i < octaves; i++) {
		unsigned int xSteps = width / octaveSize;
		unsigned int ySteps = height / octaveSize;
		unsigned int randomValueCount = (xSteps + 1) * (ySteps + 1);
		unsigned int offset;
		float *randomValues = malloc(randomValueCount * sizeof(float));
		float *xValues = malloc((width + 1) * (ySteps + 1) * sizeof(float));

		ccrSeed32(&randomizer, seed + ccnCoordinateUid(x, y) + i);

		if(makeTileable) {
			unsigned int randomValuesMinusBottom = randomValueCount - xSteps - 1;

			// Generate own
			for(j = 0; j < randomValuesMinusBottom; j++) {
				randomValues[j] = ccrGenerateFloat32(&randomizer);
			}

			// Generate bottom
			ccrSeed32(&randomizer, seed + ccnCoordinateUid(x, y + 1) + i);

			offset = (xSteps + 1) * ySteps;

			for(j = 0; j < xSteps; j++) {
				randomValues[offset + j] = ccrGenerateFloat32(&randomizer);
			}

			// Generate right
			ccrSeed32(&randomizer, seed + ccnCoordinateUid(x + 1, y) + i);

			for(j = 0; j < randomValueCount; j++) {
				unsigned int _y = j / (xSteps + 1);
				float value = ccrGenerateFloat32(&randomizer);

				if(j - _y * (xSteps + 1) == 0) randomValues[(_y + 1) * (xSteps + 1) - 1] = value;
			}

			// Generate right bottom
			ccrSeed32(&randomizer, seed + ccnCoordinateUid(x + 1, y + 1) + i);

			randomValues[randomValueCount - 1] = ccrGenerateFloat32(&randomizer);
		}
		else {
			for(j = 0; j < randomValueCount; j++) {
				randomValues[j] = ccrGenerateFloat32(&randomizer);
			}
		}

		for(j = 0; j <= ySteps ; j++) {
			for(k = 0; k <= width; k++) {
				unsigned octX = k / octaveSize;

				float factor = (float)(k - octX * octaveSize) / octaveSize;

				if(factor == 0) {
					xValues[k + j * (width + 1)] = randomValues[octX + j * (xSteps + 1)];
				}
				else {
					unsigned int index = octX + j * (xSteps + 1);

					if(interpolationMethod == CCN_INTERP_CUBIC) {
						xValues[k + j * (width + 1)] = ccTriInterpolateCubic(randomValues[octX == 0?index:index - 1], randomValues[index], randomValues[index + 1], randomValues[octX == xSteps - 1?index + 1:index + 2], factor);
					}
					else {
						xValues[k + j * (width + 1)] = ccnInterpolate(randomValues[index], randomValues[index + 1], factor, interpolationMethod);
					}
				}
			}
		}

		for(j = 0; j < size; j++) {
			unsigned int Y = j / width;
			unsigned int X = j - Y * width;

			unsigned octY = Y / octaveSize;

			float factor = (float)(Y - octY * octaveSize) / octaveSize;

			if(factor == 0) {
				(*buffer)[j] += xValues[X + octY * (width + 1)] * influence;
			}
			else {
				unsigned int index = X + octY * (width + 1);

				if(interpolationMethod == CCN_INTERP_CUBIC) {
					(*buffer)[j] += ccTriInterpolateCubic(xValues[octY == 0?index:index - width - 1], xValues[index], xValues[index + width + 1], xValues[octY == ySteps - 1?index + 1:index + ((width + 1) << 1)], factor) * influence;
				}
				else {
					(*buffer)[j] += ccnInterpolate(xValues[index], xValues[index + width + 1], factor, interpolationMethod) * influence;
				}
			}
		}
		
		free(xValues);
		free(randomValues);

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