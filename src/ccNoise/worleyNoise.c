#include <stdlib.h>
#include <assert.h>
#include <ccSort/ccSort.h>

#include "ccNoiseInternal.h"

#define _CCN_MANHATTAN_DISTANCE_FACTOR 1.414214

void ccnGenerateWorleyNoise2D(
	ccnNoise *noise,
	ccnNoiseConfiguration *configuration,
	unsigned int points,
	unsigned int n,
	unsigned int low, unsigned int high,
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
	unsigned int *pointDistances = malloc(pointListSize*sizeof(unsigned int));

	unsigned int maxManhattanDistance = (unsigned int)(high * _CCN_MANHATTAN_DISTANCE_FACTOR);

#ifdef _DEBUG
	assert(interpolationMethod != CCN_INTERP_CUBIC);
#endif

	if(configuration->tileConfiguration.tileMethod == CCN_TILE_NOT) xPeriod = yPeriod = CCN_INFINITE;

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
					pointDistances[pointId] = manhattanDistance;
				}
				else {
					pointDistances[pointId] = ccnDistance(p, pointList[j], distanceMethod);
				}

				pointId++;
			}
		}

		if(pointId > 1) ccsQuicksort(pointDistances, 0, pointId);

		if(pointId <= n || pointDistances[n] > high) {
			ccnStore(noise->values + i, configuration->storeMethod, configuration->range.high);
		}
		else if(pointDistances[n] < low) {
			ccnStore(noise->values + i, configuration->storeMethod, configuration->range.low);
		}
		else {
			ccnStore(noise->values + i, configuration->storeMethod, ccnInterpolate(configuration->range.low, configuration->range.high, (float)(pointDistances[n] - low) / (high - low), interpolationMethod));
		}
	}

	free(pointList);
	free(pointDistances);
}