#include "ccNoiseInternal.h"

void ccnGenerateWhiteNoise2D(
	ccnNoise *noise,
	ccnNoiseConfiguration *configuration)
{
	unsigned int size = noise->width * noise->height;
	unsigned int i;
	float multiplier = configuration->range.high - configuration->range.low;

	for(i = 0; i < size; ++i) {
		int Y = i / noise->width;
		int X = i - Y * noise->width;
		ccnStore(noise->values + i, configuration->storeMethod, ccrGenerateFloatCoordinate(configuration->seed, X, Y) * multiplier + configuration->range.low);
	}
}