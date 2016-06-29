#include "ccNoiseInternal.h"

void ccnGenerateWhiteNoise2D(
	ccnNoise *noise,
	const ccnNoiseConfiguration *configuration)
{
	const uint32_t size = noise->width * noise->height;
	const float multiplier = configuration->range.high - configuration->range.low;
	uint32_t i;

	for(i = 0; i < size; ++i) {
		int32_t Y = i / noise->width;
		int32_t X = i - Y * noise->width;
		ccnStore(noise->values + i, configuration->storeMethod, ccrGenerateFloatCoordinate(configuration->seed, X, Y) * multiplier + configuration->range.low);
	}
}