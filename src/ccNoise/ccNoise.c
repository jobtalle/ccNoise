#include <stdlib.h>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>
#include <ccTrigonometry/ccTrigonometry.h>
#include <ccSort/ccSort.h>

#define _CCN_CARTESIAN_INDEX(x, y, width) ((x) + (y)*(width))

#define _CCN_CENTER       0
#define _CCN_RIGHT        1
#define _CCN_RIGHT_BOTTOM 2
#define _CCN_BOTTOM       3
#define _CCN_LEFT_BOTTOM  4
#define _CCN_LEFT         5
#define _CCN_LEFT_TOP     6
#define _CCN_TOP          7
#define _CCN_RIGHT_TOP    8

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
	bool makeTileable,
	int x, int y,
	unsigned int width, unsigned int height,
	unsigned int octaves,
	unsigned int maxOctave,
	ccnInterpolationMethod interpolationMethod)
{
	unsigned int i, j;
	unsigned int size = width * height;
	unsigned int scale = maxOctave;
	
	float influence = 0.5f;

	*buffer = calloc(size, sizeof(float));

	for(i = 0; i < octaves; i++) {
		float *sourceNoise[4];

		// The dimensions for this octave
		unsigned int xSteps = width / scale;
		unsigned int ySteps = height / scale;

#define _CCN_VALUE_NOISE_SEED(X, Y) seed + ccnCoordinateUid(X, Y) + i

		ccnGenerateWhiteNoise(&sourceNoise[_CCN_CENTER], _CCN_VALUE_NOISE_SEED(x, y), xSteps, ySteps);
		ccnGenerateWhiteNoise(&sourceNoise[_CCN_RIGHT], _CCN_VALUE_NOISE_SEED(x + 1, y), xSteps, ySteps);
		ccnGenerateWhiteNoise(&sourceNoise[_CCN_RIGHT_BOTTOM], _CCN_VALUE_NOISE_SEED(x + 1, y + 1), 1, 1);
		ccnGenerateWhiteNoise(&sourceNoise[_CCN_BOTTOM], _CCN_VALUE_NOISE_SEED(x, y + 1), xSteps, 1);

#undef _CCN_VALUE_NOISE_SEED

		// Calculate Y values and add to noise
		for(j = 0; j < size; j++) {
			unsigned int Y = j / width;
			unsigned int X = j - Y * width;

			unsigned int noiseX = X / scale;
			unsigned int noiseY = Y / scale;

			float xFactor = (float)(X - noiseX * scale) / scale;
			float yFactor = (float)(Y - noiseY * scale) / scale;

			float center, right, rightBottom, bottom;

			center = sourceNoise[_CCN_CENTER][_CCN_CARTESIAN_INDEX(noiseX, noiseY, xSteps)];

			if(noiseX == xSteps - 1) {
				right = sourceNoise[_CCN_RIGHT][_CCN_CARTESIAN_INDEX(0, noiseY, xSteps)];
			}
			else {
				right = sourceNoise[_CCN_CENTER][_CCN_CARTESIAN_INDEX(noiseX + 1, noiseY, xSteps)];
			}

			if(noiseY == ySteps - 1) {
				bottom = sourceNoise[_CCN_BOTTOM][_CCN_CARTESIAN_INDEX(noiseX, 0, xSteps)];
			}
			else {
				bottom = sourceNoise[_CCN_CENTER][_CCN_CARTESIAN_INDEX(noiseX, noiseY + 1, xSteps)];
			}

			if(noiseX == xSteps - 1 && noiseY == ySteps - 1) {
				rightBottom = sourceNoise[_CCN_RIGHT_BOTTOM][0];
			}
			else {
				if(noiseX == xSteps - 1) {
					rightBottom = sourceNoise[_CCN_RIGHT][_CCN_CARTESIAN_INDEX(0, noiseY + 1, xSteps)];
				}
				else if(noiseY == ySteps - 1) {
					rightBottom = sourceNoise[_CCN_BOTTOM][_CCN_CARTESIAN_INDEX(noiseX + 1, 0, xSteps)];
				}
				else {
					rightBottom = sourceNoise[_CCN_CENTER][_CCN_CARTESIAN_INDEX(noiseX + 1, noiseY + 1, xSteps)];
				}
			}
			
			(*buffer)[j] += ccnInterpolate(
				ccnInterpolate(center, right, xFactor, interpolationMethod),
				ccnInterpolate(bottom, rightBottom, xFactor, interpolationMethod),
				yFactor, interpolationMethod) * influence;
				
		}

		// Clean up and go to next octave
		free(sourceNoise[_CCN_CENTER]);
		free(sourceNoise[_CCN_RIGHT]);
		free(sourceNoise[_CCN_RIGHT_BOTTOM]);
		free(sourceNoise[_CCN_BOTTOM]);

		scale >>= 1;
		influence *= 0.5f;
	}

	return CCN_ERROR_NONE;
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