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

static float ccnInterpolate(float a, float b, float x, ccnInterpolationMethod interpolationMethod) {
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

void ccnGenerateWorleyNoise(
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
			unsigned int manhattanDistance = abs(p.x - pointList[j].x) + abs(p.y - pointList[j].y);

			if(manhattanDistance < maxManhattanDistance) {
				switch(distanceMethod) {
				case CCN_DIST_MANHATTAN:
					pointsDistances[pointId] = manhattanDistance;
					break;
				case CCN_DIST_EUCLIDEAN:
					pointsDistances[pointId] = (int)ccTriDistance(p.x, p.y, pointList[j].x, pointList[j].y);
					break;
				case CCN_DIST_CHEBYCHEV:
					pointsDistances[pointId] = max(abs(pointList[j].x - p.x), abs(pointList[j].y - p.y));
					break;
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
}

void ccnGenerateFractalNoise(
	float **buffer,
	unsigned int seed,
	int x, int y,
	unsigned int width, unsigned int height,
	unsigned int octaves,
	unsigned int maxOctave,
	float persistence,
	ccnInterpolationMethod interpolationMethod)
{
	unsigned int size = width * height;
	unsigned int octaveSize = maxOctave;
	unsigned int i, j, k;

	float influence = 0.5f;

	ccRandomizer32 randomizer;

	octaveSize = maxOctave;

	*buffer = calloc(size, sizeof(float));

	for(i = 0; i < octaves; i++) {
		unsigned int xSteps = width / octaveSize;
		unsigned int ySteps = height / octaveSize;
		unsigned int randomValueCount = (xSteps + 1) * (ySteps + 1);
		unsigned int offset;
		float *randomValues = malloc(randomValueCount * sizeof(float));

		ccrSeed32(&randomizer, seed + ccnCoordinateUid(x, y) + i);

		// Generate block values
		for(j = 0; j < randomValueCount; j++) {
			randomValues[j] = ccrGenerateFloat32(&randomizer);
		}

		// Generate bottom
		ccrSeed32(&randomizer, seed + ccnCoordinateUid(x, y + 1) + i);

		offset = (xSteps + 1) * ySteps;

		for(j = 0; j < xSteps; j++) {
			randomValues[offset + j] = ccrGenerateFloat32(&randomizer);
		}

		// Generate right bound
		ccrSeed32(&randomizer, seed + ccnCoordinateUid(x + 1, y) + i);

		for(j = 0; j < randomValueCount; j++) {
			unsigned int _y = j / (xSteps + 1);
			float value = ccrGenerateFloat32(&randomizer);

			if(j - _y * (xSteps + 1) == 0) randomValues[(_y + 1) * (xSteps + 1) - 1] = value;
		}

		// Generate last corner value
		ccrSeed32(&randomizer, seed + ccnCoordinateUid(x + 1, y + 1) + i);

		randomValues[randomValueCount - 1] = ccrGenerateFloat32(&randomizer);

		for(j = 0; j < size; j++) {
			unsigned int Y = j / width;
			unsigned int X = j - Y * width;

			unsigned int Yoct = Y / octaveSize;
			unsigned int Xoct = (X / octaveSize) + (xSteps + 1) * Yoct;

			float xFactor = (float)(X - octaveSize * (X / octaveSize)) / octaveSize;
			float yFactor = (float)(Y - octaveSize * Yoct) / octaveSize;
			
			(*buffer)[j] += ccnInterpolate(
				ccnInterpolate(randomValues[Xoct], randomValues[Xoct + 1], xFactor, interpolationMethod),
				ccnInterpolate(randomValues[Xoct + xSteps + 1], randomValues[Xoct + xSteps + 2], xFactor, interpolationMethod),
				yFactor, interpolationMethod) * influence;
		}

		free(randomValues);

		influence *= persistence;

		octaveSize >>= 1;
		if(octaveSize == 0) break;
	}
}