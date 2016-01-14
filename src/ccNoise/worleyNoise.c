#include <stdlib.h>
#include <assert.h>

#include <ccTrigonometry/ccTrigonometry.h>

#include "ccNoiseInternal.h"

static int32_t ccnWorleyCompare(const void *a, const void *b)
{
	int32_t ia = *(int*)a;
	int32_t ib = *(int*)b;

	return ia - ib;
}

void ccnGenerateWorleyNoise2D(
	ccnNoise *noise,
	ccnNoiseConfiguration *configuration,
	uint32_t points,
	uint32_t n,
	uint32_t low, uint32_t high,
	ccnDistanceMethod distanceMethod,
	ccnInterpolationMethod interpolationMethod)
{
	uint32_t size = noise->width * noise->height;
	uint32_t i, j;
	uint32_t pointId = 0;
	uint32_t pointListSize = points * 9;
	uint32_t xPeriod = configuration->tileConfiguration.xPeriod;
	uint32_t yPeriod = configuration->tileConfiguration.yPeriod;
	ccnPoint offset;

	ccnPoint *pointList = malloc(pointListSize*sizeof(ccnPoint));
	uint32_t *pointDistances = malloc(pointListSize*sizeof(uint32_t));

	uint32_t maxManhattanDistance = (uint32_t)(high * CC_TRI_SQRT2_F);

#ifdef _DEBUG
	assert(interpolationMethod != CCN_INTERP_CUBIC);
#endif

	if(configuration->tileConfiguration.tileMethod == CCN_TILE_NOT) xPeriod = yPeriod = CCN_INFINITE;

	for(offset.x = -1; offset.x <= 1; offset.x++) {
		for(offset.y = -1; offset.y <= 1; offset.y++) {
			ccRandomizer32 randomizer;

			ccrSeed32(&randomizer, ccrGenerateUintCoordinate(configuration->seed, ccnWrapCoordinate(configuration->x + offset.x, xPeriod), ccnWrapCoordinate(configuration->y + offset.y, yPeriod)));

			for(i = 0; i < points; ++i) {
				pointList[pointId].x = (int)((ccrGenerateFloat32(&randomizer) + offset.x) * noise->width);
				pointList[pointId].y = (int)((ccrGenerateFloat32(&randomizer) + offset.y) * noise->height);
				pointId++;
			}
		}
	}

	for(i = 0; i < size; ++i) {
		ccnPoint p;
		p.y = i / noise->width;
		p.x = i - p.y * noise->width;

		pointId = 0;
		for(j = 0; j < pointListSize; j++) {
			uint32_t manhattanDistance = ccnDistance(p, pointList[j], CCN_DIST_MANHATTAN);

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

		if(pointId > 1) qsort(pointDistances, pointId, sizeof(uint32_t), ccnWorleyCompare);

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