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

	float(*interpolate)(float, float, float) = NULL;

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
			switch(interpolationMethod) {
			case CCN_INTERP_LINEAR:
				interpolate = &ccTriInterpolateLinear;
				break;
			case CCN_INTERP_QUADRATIC:
				interpolate = &ccTriInterpolateQuadratic;
				break;
			case CCN_INTERP_QUADRATIC_INVERSE:
				interpolate = &ccTriInterpolateQuadraticInverse;
				break;
			case CCN_INTERP_COSINE:
				interpolate = &ccTriInterpolateCosine;
				break;
			}

			(*buffer)[i] = interpolate(lowValue, highValue, (float)(pointsDistances[n] - low) / (high - low));
		}
	}

	free(pointList);
	free(pointsDistances);
}