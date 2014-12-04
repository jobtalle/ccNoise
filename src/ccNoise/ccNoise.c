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
	unsigned int x, y;
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
	ccnDistanceMethod distanceMethod)
{
	unsigned int size = width * height;
	unsigned int surface;
	unsigned int i, j;
	unsigned int pointId = 0;
	unsigned int pointListSize = points;

	ccnPoint *pointList = malloc(pointListSize*sizeof(ccnPoint));
	int *pointsDistances = malloc(pointListSize*sizeof(unsigned int));

	unsigned int maxManhattanDistance = high * (2 / sqrt(2));

	*buffer = malloc(sizeof(float)*size);

	for(surface = 0; surface < 1; surface++) {
		ccRandomizer32 randomizer;
		ccnPoint offset;

		offset.x = 0;
		offset.y = 0;

		ccrSeed32(&randomizer, seed + ccnCoordinateUid(x + offset.x, y + offset.y));

		for(i = 0; i < points; i++) {
			pointList[pointId].x = (int)((ccrGenerateFloat32(&randomizer) + offset.x) * width);
			pointList[pointId].y = (int)((ccrGenerateFloat32(&randomizer) + offset.y) * height);
			pointId++;
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
				}

				pointId++;
			}
		}

		if(pointId > 1) ccsQuicksort(pointsDistances, 0, pointId);
		
		if(pointId <= n) {
			(*buffer)[i] = highValue;
		}
		else if(pointsDistances[n] < low) {
			(*buffer)[i] = lowValue;
		}
		else if(pointsDistances[n] > high) {
			(*buffer)[i] = highValue;
		}
		else {
			(*buffer)[i] = lowValue + ((float)(pointsDistances[n] - low) / (high - low)) * (highValue - lowValue);
		}
	}

	free(pointList);
	free(pointsDistances);
}